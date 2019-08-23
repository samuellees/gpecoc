# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 20:18:16 2017

@author: Shone
"""


import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from utils import gol


def information_gain_ratio(Y, classes, EcocMatrix):
    """
    # take every class into consideration
    """
    num = list()
    for c in classes:
        num.append(np.sum(Y==c))

    gain_ratio = []
    for k in xrange(EcocMatrix.shape[1]):
        # make preparations
        nums_left = list()   # a list contains the number of samples for each class in left branch
        nums_right = list()   # a list contains the number of samples for each class in right branch
        for i in xrange(EcocMatrix.shape[0]):
            if EcocMatrix[i][k] == 1:
                nums_left.append(num[i])
            elif EcocMatrix[i][k] == -1:
                nums_right.append(num[i])

        # calculate global entropy
        num_left = np.float(sum(nums_left))
        num_right = np.float(sum(nums_right))
        num_total = num_left + num_right
        ps_left_global = np.array(nums_left)/num_total    # the global properties in left branch
        ps_right_global = np.array(nums_right)/num_total    # the global properties in right branch
        ent_global = -1*sum([p1*np.log2(p1) for p1 in ps_left_global]) - sum([p2*np.log2(p2) for p2 in ps_right_global])

        # calculate local entropies
        ps_left_local = np.array(nums_left)/num_left    # the local properties in left branch
        ps_right_local = np.array(nums_right)/num_right    # the local properties in right branch
        ent_left_local = -1*sum([p*np.log2(p) for p in ps_left_local])
        ent_right_local = -1*sum([p*np.log2(p) for p in ps_right_local])

        # calculate gain
        hou = (num_left*ent_left_local + num_right*ent_right_local)/num_total
        gain = ent_global - hou

        # calculate iv
        iv = -1*(num_left*np.log2(num_left/num_total) + num_right*np.log2(num_right/num_total))/num_total
        gain_ratio.append(gain/iv)
    return np.array(gain_ratio)


def information_gain(Y, classes, EcocMatrix):
    """
    # take every class into consideration
    """
    num = list()
    for c in classes:
        num.append(np.sum(Y==c))

    gain = []
    for k in xrange(EcocMatrix.shape[1]):
        # make preparations
        nums_left = list()   # a list contains the number of samples for each class in left branch
        nums_right = list()   # a list contains the number of samples for each class in right branch
        for i in xrange(EcocMatrix.shape[0]):
            if EcocMatrix[i][k] == 1:
                nums_left.append(num[i])
            elif EcocMatrix[i][k] == -1:
                nums_right.append(num[i])

        # calculate global entropy
        num_left = np.float(sum(nums_left))
        num_right = np.float(sum(nums_right))
        num_total = num_left + num_right
        ps_left_global = np.array(nums_left)/num_total    # the global properties in left branch
        ps_right_global = np.array(nums_right)/num_total    # the global properties in right branch
        ent_global = -1*sum([p*np.log2(p) for p in ps_left_global]) - sum([p*np.log2(p) for p in ps_right_global])

        # calculate local entropies
        ps_left_local = np.array(nums_left)/num_left    # the local properties in left branch
        ps_right_local = np.array(nums_right)/num_right    # the local properties in right branch
        ent_left_local = -1*sum([p*np.log2(p) for p in ps_left_local])
        ent_right_local = -1*sum([p*np.log2(p) for p in ps_right_local])

        # calculate gain
        _gain = ent_global - (num_left*ent_left_local + num_right*ent_right_local)/num_total
        gain.append(_gain)
    return np.array(gain)


def information_entropy(Y, classes, EcocMatrix):
    """
    # Calculate the Information Entropy of each column
    # At first, we need to calculate the proportion of each class
    # The result could be regard as the weight of each column
    """
    num = list()
    for c in classes:
        num.append(np.sum(Y==c))

    entropy = []
    for k in xrange(EcocMatrix.shape[1]):
        # make preparations
        nums_left = list()   # a list contains the number of samples for each class in left branch
        nums_right = list()   # a list contains the number of samples for each class in right branch
        for i in xrange(EcocMatrix.shape[0]):
            if EcocMatrix[i][k] == 1:
                nums_left.append(num[i])
            elif EcocMatrix[i][k] == -1:
                nums_right.append(num[i])

        # calculate global entropy
        num_left = np.float(sum(nums_left))
        num_right = np.float(sum(nums_right))
        num_total = num_left + num_right
        ps_left_global = np.array(nums_left)/num_total    # the global properties in left branch
        ps_right_global = np.array(nums_right)/num_total    # the global properties in right branch
        ent_global = -1*sum([p*np.log2(p) for p in ps_left_global]) - sum([p*np.log2(p) for p in ps_right_global])
        entropy.append(ent_global)
    return np.array(entropy)


def fisher_complexity(column,trainX, trainY, sel_features, thisfeature):
    """
    fisher method
    """
    feature_method_index = gol.get_val("feature_method_index")
    classes = np.array(gol.get_val("classes"))
    column = np.array(column)
    fsel_X = trainX[:, sel_features[feature_method_index[thisfeature]]]

    labels = dict()
    label_indexs = dict()
    label_indexs['1'] = np.array(np.where(column == 1))[0]
    label_indexs['-1'] = np.array(np.where(column == -1))[0]
    labels['1'] = classes[label_indexs['1']]
    bool_array = np.zeros(len(trainY), dtype=np.bool)  # generate a array being filled of 'False'
    for lb in labels['1']:
        bool_array = bool_array | np.array(trainY == lb)
    part1_X = fsel_X[bool_array]
    labels['-1'] = classes[label_indexs['-1']]
    bool_array = np.zeros(len(trainY), dtype=np.bool)  # generate a array being filled of 'False'
    for lb in labels['-1']:
        bool_array = bool_array | np.array(trainY == lb)
    part2_X = fsel_X[bool_array]

    miu1 = np.mean(np.sum(part1_X, axis=0))
    miu2 = np.mean(np.sum(part2_X, axis=0))
    sigma1 = np.var(np.sum(part1_X, axis=0))
    sigma2 = np.var(np.sum(part2_X, axis=0))

    if sigma1 + sigma2 == 0 :   return 0

    fisher_sum = (miu1 - miu2) * (miu1 - miu2) / (sigma1 + sigma2)
    return fisher_sum


def fisher_complexity_old(column,trainX, trainY, sel_features, thisfeature):
    """
    fisher method old version
    """
    feature_method_index = gol.get_val("feature_method_index")
    fea_num = gol.get_val("feature_number")
    classes = gol.get_val("classes")

    miu1 = 0; miu2 = 0; sigma1 = 0; sigma2 = 0
    d1 = np.zeros(fea_num); d2 = np.zeros(fea_num)
    d3 = np.zeros((len(classes),fea_num))    # the sum of each feature for every class
    miu3 = [0]*len(classes) # the average of d3
    sigma3 = [0]*len(classes)
    fisher11 = 0; fisher12 = 0; fisher21 = 0; fisher22 = 0
    numleft = 0; numright = 0
    num = [list(trainY).count(classes[i]) for i in xrange(len(classes))]

    # which feature selection method?
    d = trainX[:, sel_features[feature_method_index[thisfeature]]]

    # in
    for j in xrange(len(classes)):
        for m in xrange(len(trainY)):
            if (trainY[m] == classes[j]):
                d3[j,:] += d[m,:]
        if (column[j] == 1):
            numleft += num[j]
        elif (column[j] == -1):
            numright += num[j]
        miu3[j] = np.mean(d3[j])
        sigma3[j] = np.var(d3[j])


    for i in xrange(len(classes)):
        if (column[i] == 1):
            num[i] = float(num[i]) / numleft
            fisher12 += num[i] * sigma3[i]
        elif (column[i] == -1):
            num[i] = float(num[i]) / numright
            fisher22 += num[i] * sigma3[i]

    for i in xrange(len(classes)):
        for j in xrange(len(classes)):
            if (j > i):
                if (column[i] == 1 and column[j] == 1):
                    fisher11 += num[i] * num[j] * (miu3[i] - miu3[j]) * (miu3[i] - miu3[j])
                elif (column[i] == -1 and column[j] == -1):
                    fisher21 += num[i] * num[j] * (miu3[i] - miu3[j]) * (miu3[i] - miu3[j])

    # out
    for j in xrange(len(classes)):
        if column[j] == 1:
            d1[:] += d3[j,:]
        elif column[j] == -1:
            d2[:] += d3[j,:]
    miu1 = np.mean(d1)
    miu2 = np.mean(d2)
    sigma1 = np.var(d1)
    sigma2 = np.var(d2)
    if sigma1 + sigma2 == 0 :   return 0
    fisher_sum = (miu1 - miu2) * (miu1 - miu2) / (sigma1 + sigma2)

    # final
    fisherSet = []
    fisherSet.append(fisher_sum)
    if (fisher11 != 0):
        fisher1 = float(fisher11 / fisher12)
        fisherSet.append(fisher1)
    if (fisher21 != 0):
        fisher2 = float(fisher21 / fisher22)
        fisherSet.append(fisher2)

    return np.mean(fisherSet)


    

def _means_part(X,Y,labels):
    """
    # Calculate the complexity of one sungroup
    # The function can help us reduce some unnecessary code
    """
    part_X = []
    for lb in labels:
        ith_class_X = X[Y == lb]
        part_X.append(np.mean(ith_class_X, axis=0))
    euclids = euclidean_distances(part_X, part_X)
    mean_distance = np.sum(euclids)/2
    return mean_distance

def means_complexity(Train_X, Train_Y, f_m_index, f_u_list, s_features, classes, EcocMatrix):
    """
    # means, self_defined
    #  f_m_index : feature_method_index
    #   f_u_list : features_used_list
    # s_features : sel_features
    """
    # traversal every column
    classes = np.array(classes)
    Train_X = np.array(Train_X)
    Train_Y = np.array(Train_Y)

    estimators_complexity = []
    for i in xrange(EcocMatrix.shape[1]):
        labels = dict()     # show classes marked as 1,-1 and 0. eg.labels['0']=['A','C']
        label_indexs = dict()
        ith_column = EcocMatrix[:, i]
        # data for the i_th classifier
        ith_X = Train_X[:,s_features[f_m_index[f_u_list[i]]]]
        ith_Y = Train_Y
        # remove the data which are marked as '0' in the column
        label_indexs['0'] = np.array(np.where(ith_column==0))[0]
        labels['0'] = classes[label_indexs['0']]
        for lb in labels['0']:
            ith_X, ith_Y = ith_X[ith_Y!=lb], ith_Y[ith_Y!=lb]
        # divide the remaining data into two subgroups by '-1' and '1' in the column
        label_indexs['1'] = np.array(np.where(ith_column==1))[0]
        labels['1'] = classes[label_indexs['1']]
        bool_array = np.zeros(len(ith_Y),dtype = np.bool) # generate a array being filled of 'False'
        for lb in labels['1']:
            bool_array = bool_array | np.array(ith_Y==lb)
        part1_ith_X = ith_X[bool_array]
        part1_ith_Y = ith_Y[bool_array]
        
        label_indexs['-1'] = np.array(np.where(ith_column==-1))[0]
        labels['-1'] = classes[label_indexs['-1']]
        bool_array = np.zeros(len(ith_Y),dtype = np.bool) # generate a array being filled of 'False'
        for lb in labels['-1']:
            bool_array = bool_array | np.array(ith_Y==lb)
        part2_ith_X = ith_X[bool_array]
        part2_ith_Y = ith_Y[bool_array]
        # calculate the complexity between the two subgroups
        part1_ith_mean = np.mean(part1_ith_X, axis = 0)
        part2_ith_mean = np.mean(part2_ith_X, axis = 0)
        outer_ = euclidean_distances([part1_ith_mean], [part2_ith_mean])[0][0]
        
        # calculate the complexity of two subgroups respectively
        inner1_ = _means_part(part1_ith_X, part1_ith_Y, labels['1'])
        inner2_ = _means_part(part2_ith_X, part2_ith_Y, labels['-1'])
        
        # calculate complexity of current column
        estimators_complexity.append(np.mean([outer_, inner1_, inner2_]))

    # # scale to [0,1]
    # _range = [0,1]
    # _min = min(estimators_complexity)
    # _max = max(estimators_complexity)
    # return np.array([ round(((x - _min) / (1.0 * (_max - _min))) * (_range[1] - _range[0]) + _range[0],
    #               2) for x in estimators_complexity])
    return estimators_complexity
