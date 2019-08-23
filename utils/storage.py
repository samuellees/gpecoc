#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2018/1/2 12:02
@author: Eric

# 
"""
import os
import numpy as np
import cPickle as pickle
from utils import dirtools
from utils import gol

class Storager:
    def __init__(self, root_path, dataname, classifier,):
        self.root_path = root_path
        self.dataname = dataname
        self.classifier_name = classifier.__class__.__name__
        if self.classifier_name == "KNeighborsClassifier":
            self.classifier_name = "KNN"
        self.feature = None
        self.code = None

    def setfeaturecode(self, feature, code):
        self.feature = np.array(feature, dtype=int).tolist()
        self.code = np.array(code, dtype=int).tolist()

    def getfilename(self):
        str_feature = "".join(hex(i)[2:].zfill(2) for i in self.feature)
        str_code = "".join(str(i) for i in self.code)
        return str_feature + str_code + ".pkl"

    def save_estimator_train(self, estimator):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_classifier_train",
                                self.dataname, self.classifier_name, self.getfilename())
        if not os.path.exists(filepath):
            dirtools.check_folder(os.path.join(self.root_path, "_cache_", "_cache_classifier_train",
                                               self.dataname, self.classifier_name))
            with open(filepath, 'w') as f:                     # open file with write-mode
                pickle.dump(estimator, f)   # serialize and save object

    def load_estimator_train(self):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_classifier_train",
                                self.dataname, self.classifier_name, self.getfilename())
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                estimator = pickle.load(f)   # read file and build object
        else:   estimator = None
        return estimator
        # return None

    def save_estimator_test(self, estimator):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_classifier_test",
                                self.dataname, self.classifier_name, self.getfilename())
        if not os.path.exists(filepath):
            dirtools.check_folder(os.path.join(self.root_path, "_cache_", "_cache_classifier_test",
                                  self.dataname, self.classifier_name))
            with open(filepath, 'w') as f:                     # open file with write-mode
                pickle.dump(estimator, f)   # serialize and save object

    def load_estimator_test(self):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_classifier_test",
                                self.dataname, self.classifier_name, self.getfilename())
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                estimator = pickle.load(f)   # read file and build object
        else:   estimator = None
        return estimator
        # return None

    def save_prediction_valid(self, prediction):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_prediction_valid",
                                self.dataname, self.classifier_name, self.getfilename())
        if not os.path.exists(filepath):
            dirtools.check_folder(os.path.join(self.root_path, "_cache_", "_cache_prediction_valid",
                                               self.dataname, self.classifier_name))
            with open(filepath, 'w') as f:                     # open file with write-mode
                pickle.dump(prediction, f)   # serialize and save object

    def load_prediction_valid(self):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_prediction_valid",
                                self.dataname, self.classifier_name, self.getfilename())
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                prediction = pickle.load(f)   # read file and build object
        else:   prediction = None
        return prediction
        # return None

    def save_prediction_test(self, prediction):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_prediction_test",
                                self.dataname, self.classifier_name, self.getfilename())
        if not os.path.exists(filepath):
            dirtools.check_folder(os.path.join(self.root_path, "_cache_", "_cache_prediction_test",
                                               self.dataname, self.classifier_name))
            with open(filepath, 'w') as f:                     # open file with write-mode
                pickle.dump(prediction, f)   # serialize and save object

    def load_prediction_test(self):
        filepath = os.path.join(self.root_path, "_cache_", "_cache_prediction_test",
                                self.dataname, self.classifier_name, self.getfilename())
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                prediction = pickle.load(f)   # read file and build object
        else:   prediction = None
        return prediction
        # return None

