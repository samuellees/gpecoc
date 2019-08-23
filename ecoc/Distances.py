#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2018/4/10 14:07
@author: Eric

# 
"""
import numpy as np

from sklearn.utils.extmath import safe_sparse_dot
from sklearn.metrics.pairwise import check_pairwise_arrays, euclidean_distances


def get_weights(valid_output_y, ecocmatrix, valid_y):
    # prepare
    output_y_bin = np.array(valid_output_y)
    DETA = 0.01  # to avoid zero divide when calculate acccuracy
    # binarize
    for ith_output in output_y_bin:
        ith_output[ith_output > 0] = 1
        ith_output[ith_output < 0] = -1
    # compare and record the true and false number of every base_classifier for every class.
    # eg. t_num[i][j] represents the times that ith classifier truly recognized jth class.
    t_num = np.zeros((ecocmatrix.shape[1], ecocmatrix.shape[0]))
    f_num = np.zeros((ecocmatrix.shape[1], ecocmatrix.shape[0]))
    classDict = dict((j, i) for i, j in enumerate(np.unique(np.sort(valid_y))))
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
    ests_accuracy = [float(sum(t_num[i])) / (sum(t_num[i]) + sum(f_num[i]) + DETA) for i in
                     xrange(ecocmatrix.shape[1])]
    # get weights
    weights = []
    for acc in ests_accuracy:
        error = 1 - acc
        if error < 0.0000001:
            error = 0.0000001
        if error == 1:
            error = 0.9999999
        weights.append(0.5 * np.log((1 - error) / error))
    return weights


def corrected_euclidean_distances(X, Y):
    """
    X : predict values
    Y : code_book
    """
    X, Y = check_pairwise_arrays(X, Y)
    distances = safe_sparse_dot(X, Y.T, dense_output=True)
    distances *= -2
    for i in xrange(X.shape[0]):  # for each sample
        for j in xrange(Y.shape[0]):  # calculate the dist between the sample and every base classifier
            row_x = np.copy(X[i])
            row_y = np.copy(Y[j])
            row_x[row_y == 0] = 0
            distances[i][j] += np.sum(row_x * row_x)
            distances[i][j] += np.sum(row_y * row_y)
            zero_num = row_y[row_y == 0].shape[0]
            if zero_num == 0:
                distances[i][j] /= row_x.shape[0]
            else:
                distances[i][j] /= (row_x.shape[0] - zero_num)
    np.maximum(distances, 0, out=distances)
    if X is Y:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        distances.flat[::distances.shape[0] + 1] = 0.0
    return np.sqrt(distances, out=distances)


def weighting_euclidean_distances(X, Y, weights):
    """
    X : predict values
    Y : code_book
    """
    # get distance
    X, Y = check_pairwise_arrays(X, Y)
    distances = []
    for pre_vector in X:
        dists = []
        for codeword in Y:
            d1 = sum([abs(codeword[i])*weights[i]*(pre_vector[i]-codeword[i])*(pre_vector[i]-codeword[i]) for i in xrange(len(codeword))])
            dists.append(d1)
        distances.append(dists)
    distances = np.array(distances)
    distances[distances < 0] = 0
    if X is Y:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        distances.flat[::distances.shape[0] + 1] = 0.0
    return np.sqrt(distances, out=distances)


def weighting_corrected_euclidean_distances(X, Y, weights):
    """
    X : predict values
    Y : code_book
    """
    # get distance
    X, Y = check_pairwise_arrays(X, Y)
    distances = []
    for pre_vector in X:
        dists = []
        for codeword in Y:
            d1 = sum([abs(codeword[i])*weights[i]*(pre_vector[i]-codeword[i])*(pre_vector[i]-codeword[i]) for i in xrange(len(codeword))])
            d2 = sum([abs(codeword[i])*weights[i] for i in xrange(len(codeword))])
            dists.append(d1/d2)
        distances.append(dists)
    distances = np.array(distances)
    distances[distances < 0] = 0
    if X is Y:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        distances.flat[::distances.shape[0] + 1] = 0.0
    return np.sqrt(distances, out=distances)