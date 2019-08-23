#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2018/3/14 10:42
@author: Eric
Adapter for ecoc and gene bank.
# 
"""
import numpy as np
from utils import gol


# check is there exist any column being same with candidate column
def _check_duplicate(ecocmatrix, fcolumn):
    matrix_tp = ecocmatrix.transpose()
    column = np.array(fcolumn[1:])
    _bool_var = ((column == matrix_tp).all(axis=1).any() or (column*-1 == matrix_tp).all(axis=1).any())
    return _bool_var


def update_referred_times(referred_fcolumn):
    # prepare
    genebank = gol.get_val("genebank")
    # find
    for (fcolumn, est_accuracy, class_accuracies, used_frequence) in genebank.genes:
        if (referred_fcolumn == fcolumn).all() or (referred_fcolumn[0] == fcolumn[0] and referred_fcolumn[1:]*-1 == fcolumn[1:]).all():
            # update
            used_frequence[0] += 1
            break


def get_gene_from_bank(features, ecocmatrix, confusion_matrix, classes):
    # prepare
    genebank = gol.get_val("genebank")
    DETA = 0.01  # to avoid zero divide when calculate acccuracy
    # calculate the errors
    errors = list()
    for i in xrange(len(classes)):
        errors.append(1 - float(confusion_matrix[i, i]) / np.sum(confusion_matrix[i, :]))
    # select classes that have higher error than average.
    avg_errors = np.mean(errors)
    hard_classes = np.where(errors > avg_errors)[0]
    if len(hard_classes) == 0: return # it means accucacy is 1
    # calculate score of every gene, and sort by decs.
    score_tuple = []
    for (fcolumn, est_accuracy, class_accuracies, used_frequence) in genebank.genes:
        scores = [(errors[i] - (1 - class_accuracies[i])) * abs(fcolumn[i+1]) for i in hard_classes]
        score = sum(scores)/(sum([abs(fcolumn[i+1]) + DETA for i in hard_classes]))
        score_tuple.append((fcolumn, score))
    score_tuple = sorted(score_tuple, key=lambda s: s[1], reverse=True)
    # select most suitable column
    candidate = None
    for (fcolumn, score) in score_tuple:
        if not _check_duplicate(ecocmatrix, fcolumn):
            candidate = fcolumn
            break
    return candidate


def add_gene_from_matrix(features, ecocmatrix, output_y, valid_y, classes):
    # prepare
    ADD_PERCENT = 0.3   # the percent of columns to be saved.
    DETA = 0.01  # to avoid zero divide when calculate acccuracy
    genebank = gol.get_val("genebank")
    output_y_bin = np.array(output_y)
    # binarize
    for ith_output in output_y_bin:
        ith_output[ith_output > 0] = 1
        ith_output[ith_output < 0] = -1
    # compare and record the true and false number of every base_classifier for every class.
    # eg. t_num[i][j] represents the times that ith classifier truly recognized jth class.
    t_num = np.zeros((ecocmatrix.shape[1], ecocmatrix.shape[0]))
    f_num = np.zeros((ecocmatrix.shape[1], ecocmatrix.shape[0]))
    classDict = dict((j, i) for i, j in enumerate(classes))
    for i in xrange(len(output_y_bin)):
        ith_output = output_y_bin[i]
        ith_codeword = ecocmatrix[classDict[valid_y[i]]]
        for j in xrange(len(ith_codeword)):
            if ith_codeword[j] == 0: continue
            if ith_codeword[j] == ith_output[j]:
                t_num[j][classDict[valid_y[i]]] += 1
            else:
                f_num[j][classDict[valid_y[i]]] += 1
    # calculate the accuracy of every base_classifier.
    ests_accuracy = [float(sum(t_num[i]))/(sum(t_num[i])+sum(f_num[i])+DETA) for i in xrange(ecocmatrix.shape[1])]
    # calculate accuracy of every base_classifier for every class.
    est_class_accuracies = []
    for i in xrange(ecocmatrix.shape[1]):
        est_class_accuracies.append([(i, j, float(t_num[i][j])/(t_num[i][j]+f_num[i][j]+DETA)) for j in xrange(ecocmatrix.shape[0])])
    # save n column randomly to genebank
    n = np.ceil(ecocmatrix.shape[1] * ADD_PERCENT)
    for i in xrange(int(n)):
        import random
        est_index = random.randint(0, ecocmatrix.shape[1]-1)
        est_accuracy = ests_accuracy[est_index]
        feature = features[est_index]
        col = ecocmatrix[:, est_index]
        fcolumn = np.hstack((feature, col))
        class_accuracies = [accuracy for (i_est, j_cls, accuracy) in est_class_accuracies[est_index]]
        genebank.addgene(fcolumn, est_accuracy, class_accuracies, [0])







