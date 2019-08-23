# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 11:00:23 2017

@author: Eric
"""
import numpy as np
import gp.TreeMatrixConvertor as TMConvertor
from utils import gol
from scipy.spatial import distance
from ConnectClassifier import ConnectClassifier as CC
from sklearn.metrics.pairwise import euclidean_distances
from DataComplexity import means_complexity
from DataComplexity import information_gain
from DataComplexity import information_gain_ratio
from DataComplexity import information_entropy

def eval_func_test(chromosome):
    """
    # calculate fscore
    """
    print type(chromosome)
    None


def eval_func_fscore(chromosome):
    """
    # calculate fscore
    """
    EcocMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(chromosome)
    Train_X = gol.get_val("Train_X")
    Train_Y = gol.get_val("Train_Y")
    validation_X = gol.get_val("validation_X")
    validation_Y = gol.get_val("validation_Y")
    sel_features = gol.get_val("sel_features")

    cc = CC(features_used_list, sel_features, EcocMatrix)
    fscore, accuracy, infos_evaluations = cc.TrainAndTest(Train_X, Train_Y, validation_X, validation_Y)
    chromosome.infos_evaluation = infos_evaluations
    return fscore, accuracy


def eval_func_eucdist(chromosome):
    """
    # calculate avg_euclidean_dist of a individual
    """
    EcocMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(chromosome)
    classes = gol.get_val("classes")
    num_class = len(classes)
    num_cols = EcocMatrix.shape[1]
    _dist = euclidean_distances(EcocMatrix, EcocMatrix)/np.sqrt(num_cols)
    dist = np.sum(_dist)/2/(num_class*(num_class-1))
    return dist


def eval_func_entropy(chromosome):
    """
    # Calculate the complexity named "means"
    # The data is all training set
    """
    Train_Y = gol.get_val("Train_Y")
    classes = gol.get_val("classes")
    EcocMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(chromosome)
    entropy = information_entropy(Train_Y, classes, EcocMatrix)
    return np.mean(entropy)


def eval_func_information_gain(chromosome):
    """
    # Calculate the information gain
    # The data is all training set
    """
    Train_Y = gol.get_val("Train_Y")
    classes = gol.get_val("classes")
    EcocMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(chromosome)
    infor_gain = information_gain(Train_Y, classes, EcocMatrix)
    return np.mean(infor_gain)


def eval_func_hamdist(chromosome):
    """
    # calculate hamdist of a individual
    """
    EcocMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(chromosome)
    classes = gol.get_val("classes")
    dist = 0
    for i in xrange(len(EcocMatrix)):
        for j in xrange(i+1, len(EcocMatrix)):
            dist += distance.hamming(EcocMatrix[i], EcocMatrix[j])
    num = len(classes)*(len(classes)-1)/2
    dist /= num
    return dist

