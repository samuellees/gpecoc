#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 20:47:46 2017

@author: hanrui
"""

from utils import gol
from OutputCodeClassifier import OutputCodeClassifier

import numpy as np

from sklearn.naive_bayes import GaussianNB,BernoulliNB
from sklearn.svm import LinearSVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier


class ConnectClassifier:

    def __init__(self, features_used_list, sel_features, codeMatrix):
        self.sel_features = sel_features
        self.code_Matrix = np.array(codeMatrix)
        self.features_used_list = features_used_list
        self.estimator = LinearSVC(random_state=0)
        # self.estimator = DecisionTreeClassifier(random_state=0)
        # self.estimator = GaussianNB()
        # self.estimator = BernoulliNB()
        # self.estimator = KNeighborsClassifier(n_neighbors=gol.get_val("n_neighbors"))
        self.oc = OutputCodeClassifier(self.estimator, random_state=0)


    def TrainAndTest(self, train_X, train_Y, validation_X, validation_Y):
        self.oc.fit(train_X, train_Y, self.features_used_list, self.sel_features, self.code_Matrix)
        score, accuracy, infos_evaluations = self.oc.predict(self.features_used_list, self.sel_features,
                                        train_X, train_Y, validation_X, validation_Y)
        return score, accuracy, infos_evaluations


    def FinalTrainAndTest(self, train_X, train_Y, validation_X, validation_Y, test_X, test_Y):
        self.oc.fit(train_X, train_Y, self.features_used_list,  self.sel_features, self.code_Matrix)
        final_score, final_accuracy, infos_evaluations = self.oc.predictFinal(self.features_used_list, self.sel_features,
                                                        train_X, train_Y, validation_X, validation_Y, test_X, test_Y)
        return final_score, final_accuracy, infos_evaluations


    def TrainAndTest_withoutlocalimp(self, train_X, train_Y, validation_X, validation_Y):
        self.oc.fit(train_X, train_Y, self.features_used_list, self.sel_features, self.code_Matrix)
        score, accuracy = self.oc.predict_withoutlocalimprovement(self.features_used_list, self.sel_features,
                                        train_X, train_Y, validation_X, validation_Y)
        return score, accuracy


    def FinalTrainAndTest_withoutlocalimp(self, train_X, train_Y, validation_X, validation_Y, test_X, test_Y):
        self.oc.fit(train_X, train_Y, self.features_used_list,  self.sel_features, self.code_Matrix)
        final_score, final_accuracy = self.oc.predictFinal_withoutlocalimprovement(self.features_used_list, self.sel_features,
                                                        train_X, train_Y, validation_X, validation_Y, test_X, test_Y)
        return final_score, final_accuracy
