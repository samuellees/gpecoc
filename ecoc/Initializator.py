# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:19:59 2017

@author: Shone
"""
import math

import numpy as np
import Configurations as Configs

from preprocess import DataLoader, FeatureSelection
from utils import gol
from utils import genebank as gbank

"""
# init all thing
"""
def init_all():
    # init_gol()
    init_config()
    init_dataset()
    init_classes()
    init_feature()
    init_maxColumn_and_maxDeap()
    init_genebank()


"""
# init module 'gol'
"""
def init_gol():
    gol._init()

"""
# load 'config' in gol
"""
def init_config():
    dataName = gol.get_val("dataName")
    gol.set_val("n_jobs", Configs.n_jobs)
    gol.set_val("version", Configs.version)
    gol.set_val("testFile", "data/" + dataName + "_test.data")
    gol.set_val("trainFile", "data/" + dataName + "_train.data")
    gol.set_val("validationFile", "data/" + dataName + "_validation.data")
    gol.set_val("root_path", Configs.root_path)
    gol.set_val("growMethod", Configs.growMethod)
    gol.set_val("freq_stats", Configs.freq_stats)
    gol.set_val("generations", Configs.generations)
    gol.set_val("n_neighbors", Configs.n_neighbors)
    gol.set_val("mutationRate", Configs.mutationRate)
    gol.set_val("crossoverRate", Configs.crossoverRate)
    gol.set_val("populationSize", Configs.populationSize)

"""
# process and load dataset in gol
"""
def init_dataset():
    trainFile = gol.get_val("trainFile")
    testFile = gol.get_val("testFile")
    validationFile = gol.get_val("validationFile")
    Train_X, Train_Y, validation_X, validation_Y, Test_X , Test_Y = DataLoader.loadDataset(trainFile, validationFile, testFile)
    class_num = len(np.unique(Train_Y))
    length = len(Train_Y) + len(validation_Y) + len(Test_Y)
    gol.set_val("Train_X", Train_X)
    gol.set_val("Train_Y", Train_Y)
    gol.set_val("validation_X", validation_X)
    gol.set_val("validation_Y", validation_Y)
    gol.set_val("Test_X", Test_X)
    gol.set_val("Test_Y", Test_Y)

"""
# load classes in gol
"""
def init_classes():
    trainFile = gol.get_val("trainFile")
    validationFile = gol.get_val("validationFile")
    testFile = gol.get_val("testFile")
    classes = DataLoader.loadClasses(trainFile, validationFile, testFile)   # a list
    gol.set_val("classes", classes)

"""
# load all thing about feature in gol
"""
def init_feature():
    # features
    Train_X = gol.get_val("Train_X")
    Train_Y = gol.get_val("Train_Y")
    trainFile = gol.get_val("trainFile")
    # the num of feature to be selected is half of the feature numbers
    feature_number = Train_X.shape[1]/2+1 if Train_X.shape[1]/2+1<75 else 75
    fea = FeatureSelection.select_features(trainFile, Train_X, Train_Y, feature_number)
    sel_features_backup = fea[0]
    feature_F1 = fea[1]
    feature_F2 = fea[2]
    feature_F3 = fea[3]
    feature_F4 = fea[4]
    sel_features = []
    sel_features.append(feature_F1)
    sel_features.append(feature_F2)
    sel_features.append(feature_F3)
    sel_features.append(feature_F4)
    # feature_method
    # feature_method_index
    feature_method_index = dict((c, i) for i, c in enumerate(FeatureSelection.feature_method))
    gol.set_val("sel_features",sel_features)
    gol.set_val("feature_number", feature_number)
    gol.set_val("feature_method_index",feature_method_index)
    gol.set_val("feature_method",FeatureSelection.feature_method)

"""
# init maxcolumn and maxdeap
# Designed following the theory "2N_c"
"""
def init_maxColumn_and_maxDeap():
    n_classes = float(len(gol.get_val("classes")))
    maxColumn = n_classes * 2
    maxDeap = np.ceil(math.log(maxColumn,2))
    gol.set_val("maxColumn",int(maxColumn))
    gol.set_val("maxDeap",int(maxDeap))

"""
# init genebank to save good column
"""
def init_genebank():
    genebank = gbank.GeneBank()
    gol.set_val("genebank", genebank)

