#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# float /

from __future__ import division
from sklearn.metrics import confusion_matrix
from sklearn.utils.validation import check_is_fitted
from sklearn.base import MetaEstimatorMixin
from sklearn.metrics import precision_recall_fscore_support
from sklearn.base import BaseEstimator, ClassifierMixin, clone

import numpy as np
import warnings
import copy, random
import utils.LegalityCheckers as LC

from ecoc.Distances import corrected_euclidean_distances, \
    weighting_corrected_euclidean_distances, \
    weighting_euclidean_distances, get_weights, euclidean_distances
from ecoc.Utils import get_gene_from_bank, add_gene_from_matrix, update_referred_times
from utils import gol, delog
from utils.storage import Storager



def _fit_binary(estimator, X, y, classes=None):
    """Fit a single binary estimator."""
    unique_y = np.unique(y)
    if len(unique_y) == 1:
        if classes is not None:
            if y[0] == -1:
                c = 0
            else:
                c = y[0]
            warnings.warn("Label %s is present in all training examples." %
                          str(classes[c]))
        estimator = _ConstantPredictor().fit(X, unique_y)
    else:
        estimator = clone(estimator)
        estimator.fit(X, y)
    return estimator


def corrected_fit_binary(estimator, X, y, classes=None):
    """Fit a single binary estimator."""
    X, y = X[y != 0], y[y != 0]
    y[y == -1] = 0
    return _fit_binary(estimator, X, y)


def corrected_predict_binary(estimator, X):
    """Make predictions using a single binary estimator."""
    score = _predict_binary(estimator, X)
    score *= 2
    score -= 1
    return np.array(score, dtype=np.float)


def _partial_fit_binary(estimator, X, y):
    """Partially fit a single binary estimator."""
    estimator.partial_fit(X, y, np.array((0, 1)))
    return estimator

def _predict_binary(estimator, X):
    """Make predictions using a single binary estimator."""
    if getattr(estimator, "_estimator_type", None) == "regressor":
        return estimator.predict(X)
    try:
        score = np.ravel(estimator.decision_function(X))
    except (AttributeError, NotImplementedError):
        # probabilities of the positive class
        score = estimator.predict_proba(X)[:, 1]

    return score


def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if (not hasattr(estimator, "decision_function") and
            not hasattr(estimator, "predict_proba")):
        raise ValueError("The base estimator should implement "
                         "decision_function or predict_proba!")


def _sigmoid_normalize(X):  # sigmoid
    X = (X+1)/2
    return (1 / (1 + np.exp(-X)))*2-1


class _ConstantPredictor(BaseEstimator):
    def fit(self, X, y):
        self.y_ = y
        return self

    def predict(self, X):
        check_is_fitted(self, 'y_')
        return np.repeat(self.y_, X.shape[0])

    def decision_function(self, X):
        check_is_fitted(self, 'y_')
        return np.repeat(self.y_, X.shape[0])

    def predict_proba(self, X):
        check_is_fitted(self, 'y_')
        return np.repeat([np.hstack([1 - self.y_, self.y_])],
                         X.shape[0], axis=0)


class OutputCodeClassifier(BaseEstimator, ClassifierMixin, MetaEstimatorMixin):
    def __init__(self, estimator, random_state=None):
        self.classes = gol.get_val("classes")
        self.classes_ = None
        self.code_book_ = None
        self.estimator = estimator
        self.estimators_ = None
        self.estimator_type = None
        self.featuresNames = gol.get_val("feature_method")
        self.infos_evaluations = []
        self.n_jobs = gol.get_val("n_jobs")
        self.random_state = random_state
        self.storager = Storager(gol.get_val("root_path"), gol.get_val("dataName"), self.estimator)
        self.trainX = None
        self.weights = None

    def fit(self, X, y, features_used_list, sel_features, code_book):

        _check_estimator(self.estimator)
        if hasattr(self.estimator, "decision_function"):
            self.estimator_type = 'decision_function'  # output = [-Nan,Nan]
        else:
            self.estimator_type = 'predict_proba'  #  output = [0, 1]

        self.classes_ = np.unique(np.sort(y))
        self.code_book_ = code_book
        self.trainX = X
        feature_method_index = gol.get_val("feature_method_index")
        classes_index = dict((c, i) for i, c in enumerate(self.classes_))
        extend_ecocmatrix = np.array([self.code_book_[classes_index[y[i]]]
                      for i in range(X.shape[0])], dtype=np.int)
        # try to restore estimators from cache
        self.estimators_ = list()
        for i in range(code_book.shape[1]):
            _column = self.code_book_[:, i]
            _features = feature_method_index[features_used_list[i]]
            self.storager.setfeaturecode(sel_features[_features], _column)
            est = self.storager.load_estimator_train()
            if est is None:
                # need training
                est = corrected_fit_binary(self.estimator, X[:, sel_features[_features]], extend_ecocmatrix[:, i])
                self.storager.save_estimator_train(est)
            self.estimators_.append(est)
        return self

    def fit_one(self, train_x, train_y, sel_features, fcolumn):
        # prepare
        feature = fcolumn[0]
        column = fcolumn[1:]
        classes_index = dict((c, i) for i, c in enumerate(self.classes_))
        # try to store from cache
        self.storager.setfeaturecode(sel_features[feature], column)
        new_estimator = self.storager.load_estimator_train()
        if new_estimator is None:
            # training
            extend_column = np.array([column[classes_index[train_y[i]]]
                                for i in xrange(train_x.shape[0])], dtype=np.int)
            new_estimator = corrected_fit_binary(self.estimator,
                                train_x[:, sel_features[feature]], extend_column)
            self.storager.save_estimator_train(new_estimator)
        return new_estimator

    def predict(self, features_used_list, sel_features, Train_X, Train_Y, Valid_X, valid_Y):
        self.feature_name = features_used_list
        self.trainX = Train_X
        self.trainY = Train_Y
        self.sel_features = sel_features
        # prepare
        feature_method_index = gol.get_val("feature_method_index")
        check_is_fitted(self, 'estimators_')
        # try to restore output from cache
        output_y = []
        for i in xrange(len(self.estimators_)):
            _column = self.code_book_[:, i]
            _features = feature_method_index[features_used_list[i]]
            self.storager.setfeaturecode(sel_features[_features], _column)
            pre = self.storager.load_prediction_valid()
            if pre is None:
                # need predicting
                pre = corrected_predict_binary(self.estimators_[i], Valid_X[:, sel_features[_features]])
                self.storager.save_prediction_valid(pre)
            output_y.append(pre)
        output_y = np.array(output_y).T
        if self.estimator_type == 'decision_function':
            output_y = _sigmoid_normalize(output_y)  
        # get score and confusion matrix
        pred = self.get_distances(output_y, self.code_book_).argmin(axis=1)
        score, accuracy = self.calculateFScore(self.classes_[pred], valid_Y)
        self.conMatrix = confusion_matrix(valid_Y, self.classes_[pred])
        # infos about features
        _message = "Performance without local improvement:"
        self.matrix_tracer(_message, score, accuracy)
        test_X = gol.get_val("Test_X")
        test_Y = gol.get_val("Test_Y")
        t_score, t_acc = self.predictFinal_withoutlocalimprovement(features_used_list, self.sel_features,
                                                                   Train_X, Train_Y, Valid_X, valid_Y, test_X, test_Y)
        self.infos_evaluations.insert(len(self.infos_evaluations), "test-f-score:"+str(t_score))
        self.infos_evaluations.insert(len(self.infos_evaluations), "test-accuracy:"+str(t_acc))
        self.infos_evaluations.insert(len(self.infos_evaluations), self.conMatrix)

        # adding column
        temp_features_name = dict((c, i) for (i, c) in feature_method_index.items())
        add_counter = 0
        while True:
            features_digit = [feature_method_index[features_used_list[i]] for i in xrange(self.code_book_.shape[1])]
            add_fcol = get_gene_from_bank(features_digit, self.code_book_, self.conMatrix, self.classes_)
            if add_fcol is not None:
                # prepare new_ecocmatrix and output_y
                new_ecocmatrix = np.hstack((self.code_book_, np.array([add_fcol[1:]]).transpose()))
                # reconstruct the output_y because they can not be sigmoid respectively
                new_y = []
                for i in xrange(len(self.estimators_)):
                    _column = self.code_book_[:, i]
                    _features = feature_method_index[features_used_list[i]]
                    self.storager.setfeaturecode(sel_features[_features], _column)
                    pre = self.storager.load_prediction_valid()
                    if pre is None:
                        # need predicting
                        pre = corrected_predict_binary(self.estimators_[i], Valid_X[:, sel_features[_features]])
                        self.storager.save_prediction_valid(pre)
                    new_y.append(pre)
                # add new output of new column, need training
                new_estimator = self.fit_one(Train_X, Train_Y, sel_features, add_fcol)
                add_fcol_y = corrected_predict_binary(new_estimator, Valid_X[:, sel_features[add_fcol[0]]])
                new_y.append(add_fcol_y)
                new_y = np.array(new_y).T
                if self.estimator_type == 'decision_function':
                    new_y = _sigmoid_normalize(new_y)
                # calculate new accuracy and compare
                new_pred = self.get_distances(new_y, new_ecocmatrix).argmin(axis=1)
                new_score, new_accuracy = self.calculateFScore(self.classes_[new_pred], valid_Y)
                # update if there be any improvement
                if new_accuracy > accuracy:
                # if new_score >= score:
                    output_y = new_y
                    pred = new_pred
                    score = new_score
                    accuracy = new_accuracy
                    self.code_book_ = new_ecocmatrix
                    self.feature_name = np.hstack((features_used_list, temp_features_name[add_fcol[0]]))
                    features_used_list = self.feature_name
                    self.estimators_.insert(len(self.estimators_), new_estimator)
                    self.conMatrix = confusion_matrix(valid_Y, self.classes_[pred])
                    # update used frequency
                    update_referred_times(add_fcol)
                    # update matrix tracer
                    add_counter += 1
                    _message = str(add_counter) + "Add one column:"
                    self.matrix_tracer(_message, score, accuracy)
                    test_X = gol.get_val("Test_X")
                    test_Y = gol.get_val("Test_Y")
                    t_score, t_acc = self.predictFinal_withoutlocalimprovement(features_used_list, self.sel_features,
                                                                               Train_X, Train_Y, Valid_X, valid_Y, test_X, test_Y)
                    self.infos_evaluations.insert(len(self.infos_evaluations), "test-f-score:"+str(t_score))
                    self.infos_evaluations.insert(len(self.infos_evaluations), "test-accuracy:"+str(t_acc))
                    self.infos_evaluations.insert(len(self.infos_evaluations), self.conMatrix)
                else:
                    # no improvement, stop.
                    break
            else:
                # genebank is empty, or no suitable column, stop.
                break

        # save good column into genebank
        features_digit = [feature_method_index[features_used_list[i]] for i in xrange(self.code_book_.shape[1])]
        add_gene_from_matrix(features_digit, self.code_book_, output_y, valid_Y, self.classes_)
        return score, accuracy, self.infos_evaluations

    def predictFinal(self, features_used_list, sel_features, train_X, train_Y, valid_X, valid_Y, test_X, test_Y):
        self.feature_name = features_used_list
        self.trainY = train_Y
        self.sel_features = sel_features
        # prepare
        check_is_fitted(self, 'estimators_')
        feature_method_index = gol.get_val("feature_method_index")
        # try to restore output from cache
        output_y = []
        for i in xrange(len(self.estimators_)):
            _column = self.code_book_[:, i]
            _features = feature_method_index[features_used_list[i]]
            self.storager.setfeaturecode(sel_features[_features], _column)
            pre = self.storager.load_prediction_valid()
            if pre is None:
                pre = corrected_predict_binary(self.estimators_[i], valid_X[:, sel_features[_features]])
                self.storager.save_prediction_valid(pre)
            output_y.append(pre)
        output_y = np.array(output_y).T
        if self.estimator_type == 'decision_function':
            output_y = _sigmoid_normalize(output_y)  
        # get score and confusion matrix
        pred = self.get_distances(output_y, self.code_book_).argmin(axis=1)
        score, accuracy = self.calculateFScore(self.classes_[pred], valid_Y)
        self.conMatrix = confusion_matrix(valid_Y, self.classes_[pred])
        # log
        _message = "Performance without local improvement:"
        self.matrix_tracer(_message, score, accuracy)
        t_score, t_acc = self.predictFinal_withoutlocalimprovement(features_used_list, self.sel_features,
                                                                   train_X, train_Y, valid_X, valid_Y, test_X, test_Y)
        self.infos_evaluations.insert(len(self.infos_evaluations), "test-f-score:"+str(t_score))
        self.infos_evaluations.insert(len(self.infos_evaluations), "test-accuracy:"+str(t_acc))
        self.infos_evaluations.insert(len(self.infos_evaluations), self.conMatrix)

        # adding column
        temp_features_name = dict((c, i) for i, c in feature_method_index.items())
        add_counter = 0
        while True:
            features_digit = [feature_method_index[features_used_list[i]] for i in xrange(self.code_book_.shape[1])]
            add_fcol = get_gene_from_bank(features_digit, self.code_book_, self.conMatrix, self.classes_)
            if add_fcol is not None:
                # prepare new_ecocmatrix and output_y
                new_ecocmatrix = np.hstack((self.code_book_, np.array([add_fcol[1:]]).transpose()))
                # reconstruct the output_y because they can not be sigmoid respectively
                new_y = []
                for i in xrange(len(self.estimators_)):
                    _column = self.code_book_[:, i]
                    _features = feature_method_index[features_used_list[i]]
                    self.storager.setfeaturecode(sel_features[_features], _column)
                    pre = self.storager.load_prediction_valid()
                    if pre is None:
                        pre = corrected_predict_binary(self.estimators_[i], valid_X[:, sel_features[_features]])
                        self.storager.save_prediction_valid(pre)
                    new_y.append(pre)
                # add new output of new column, need training
                new_estimator = self.fit_one(train_X, train_Y, sel_features, add_fcol)
                add_fcol_y = corrected_predict_binary(new_estimator, valid_X[:, sel_features[add_fcol[0]]])
                new_y.append(add_fcol_y)
                new_y = np.array(new_y).T
                if self.estimator_type == 'decision_function':
                    new_y = _sigmoid_normalize(new_y)
                # calculate new accuracy and compare
                new_pred = self.get_distances(new_y, new_ecocmatrix).argmin(axis=1)
                new_score, new_accuracy = self.calculateFScore(self.classes_[new_pred], valid_Y)
                # update if there be any improvement
                if new_accuracy > accuracy:
                # if new_score >= score:
                    output_y = new_y
                    pred = new_pred
                    score = new_score
                    accuracy = new_accuracy
                    self.code_book_ = new_ecocmatrix
                    self.feature_name = np.hstack((features_used_list, temp_features_name[add_fcol[0]]))
                    features_used_list = self.feature_name
                    self.estimators_.insert(len(self.estimators_), new_estimator)
                    self.conMatrix = confusion_matrix(valid_Y, self.classes_[pred])
                    # update used frequency
                    update_referred_times(add_fcol)
                    # update matrix tracer
                    add_counter += 1
                    _message = str(add_counter) + "Add one column:"
                    self.matrix_tracer(_message, score, accuracy)
                    t_score, t_acc = self.predictFinal_withoutlocalimprovement(features_used_list, self.sel_features,
                                                                               train_X, train_Y, valid_X, valid_Y, test_X, test_Y)
                    self.infos_evaluations.insert(len(self.infos_evaluations), "test-f-score:"+str(t_score))
                    self.infos_evaluations.insert(len(self.infos_evaluations), "test-accuracy:"+str(t_acc))
                    self.infos_evaluations.insert(len(self.infos_evaluations), self.conMatrix)
                else:
                    # no improvement, stop.
                    break
            else:
                # genebank is empty, or no suitable column, stop.
                break
        ########
        # TEST #
        ########
        # retraining because different training set, try to restore estimators from cache.
        classes_index = dict((c, i) for i, c in enumerate(self.classes_))
        final_train_x = np.vstack((train_X, valid_X))
        final_train_y = np.hstack((train_Y, valid_Y))
        self.estimators_ = list()
        for i in range(self.code_book_.shape[1]):
            _column = self.code_book_[:, i]
            _features = feature_method_index[features_used_list[i]]
            self.storager.setfeaturecode(sel_features[_features], _column)
            est = self.storager.load_estimator_test()
            if est is None:
                # need training
                extend_column = np.array([_column[classes_index[final_train_y[i]]]
                                    for i in xrange(final_train_x.shape[0])], dtype=np.int)
                est = corrected_fit_binary(self.estimator, final_train_x[:, sel_features[_features]], extend_column)
                self.storager.save_estimator_test(est)
            self.estimators_.append(est)
        # predicting because different training set, try to restore output from cache.
        output_y = []
        for i in xrange(len(self.estimators_)):
            _column = self.code_book_[:, i]
            _features = feature_method_index[features_used_list[i]]
            self.storager.setfeaturecode(sel_features[_features], _column)
            pre = self.storager.load_prediction_test()
            if pre is None:
                pre = corrected_predict_binary(self.estimators_[i], test_X[:, sel_features[_features]])
                self.storager.save_prediction_test(pre)
            output_y.append(pre)
        output_y = np.array(output_y).T
        if self.estimator_type == 'decision_function':
            output_y = _sigmoid_normalize(output_y)  
        # get score
        pred = self.get_distances(output_y, self.code_book_, Weighted=True).argmin(axis=1)
        score, accuracy = self.calculateFScore(self.classes_[pred], test_Y)
        return score, accuracy, self.infos_evaluations

    def calculateFScore(self, pred, validation_labels):
        # calculate f-score
        accuracy = (np.float)(pred[pred == validation_labels].shape[0]) / (np.float)(pred.shape[0])
        scoreList = precision_recall_fscore_support(validation_labels, pred)[2]
        score = np.average(scoreList)
        return score, accuracy


    def get_distances(self, output_y, code_book_, Weighted=False):
        # need weighted
        if not Weighted:
            valid_y = gol.get_val("validation_Y")
            self.weights = get_weights(output_y, code_book_, valid_y)
        return weighting_corrected_euclidean_distances(output_y, code_book_, self.weights)


    def matrix_tracer(self, _message, score, accuracy):
        self.infos_evaluations.insert(len(self.infos_evaluations), _message)
        self.infos_evaluations.insert(len(self.infos_evaluations), self.feature_name)
        self.infos_evaluations.insert(len(self.infos_evaluations), self.code_book_)
        self.infos_evaluations.insert(len(self.infos_evaluations), 'Confusion Matrix:')
        self.infos_evaluations.insert(len(self.infos_evaluations), self.conMatrix)
        self.infos_evaluations.insert(len(self.infos_evaluations), "f-score:"+str(score))
        self.infos_evaluations.insert(len(self.infos_evaluations), "accuracy:"+str(accuracy)+'\n')

    def predict_withoutlocalimprovement(self, features_used_list, sel_features, Train_X, Train_Y, Valid_X, valid_Y):
        self.feature_name = features_used_list
        self.trainX = Train_X
        self.trainY = Train_Y
        self.sel_features = sel_features

        check_is_fitted(self, 'estimators_')
        feature_method_index = gol.get_val("feature_method_index")

        Y = []
        for i in xrange(len(self.estimators_)):
            self.storager.setfeaturecode(sel_features[feature_method_index[features_used_list[i]]], self.code_book_[:, i])
            pre = self.storager.load_prediction_valid()
            if pre is None:
                pre = corrected_predict_binary(self.estimators_[i],
                                               Valid_X[:, sel_features[feature_method_index[features_used_list[i]]]])
                self.storager.save_prediction_valid(pre)
            Y.append(pre)
        Y = np.array(Y).T

        if self.estimator_type == 'decision_function':
            Y = _sigmoid_normalize(Y)  
        pred = self.get_distances(Y, self.code_book_).argmin(axis=1)
        self.conMatrix = confusion_matrix(valid_Y, self.classes_[pred])

        score, accuracy = self.calculateFScore(self.classes_[pred], valid_Y)
        return score, accuracy

    def predictFinal_withoutlocalimprovement(self, features_used_list, sel_features, train_X, train_Y, valid_X, valid_Y, test_X, test_Y):

        self.feature_name = features_used_list
        self.trainX = np.append(train_X, valid_X, axis=0)
        self.trainY = np.append(train_Y, valid_Y)
        self.sel_features = sel_features

        feature_method_index = gol.get_val("feature_method_index")

        self.fit(self.trainX, self.trainY, features_used_list, sel_features, self.code_book_)
        check_is_fitted(self, 'estimators_')

        Y = []
        for i in xrange(len(self.estimators_)):
            pre = corrected_predict_binary(self.estimators_[i],
                                               test_X[:, sel_features[feature_method_index[features_used_list[i]]]])
            Y.append(pre)
        Y = np.array(Y).T
        if self.estimator_type == 'decision_function':
            Y = _sigmoid_normalize(Y) 
        pred = self.get_distances(Y, self.code_book_, Weighted=True).argmin(axis=1)
        self.conMatrix = confusion_matrix(test_Y, self.classes_[pred])
        score, accuracy = self.calculateFScore(self.classes_[pred], test_Y)
        return score, accuracy
