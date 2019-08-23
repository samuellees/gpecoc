#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:42:04 2017

@author: root
"""

import numpy as np
from sklearn.model_selection import StratifiedKFold


'''
# split given data into several kolds
'''
class KFold:
     
    def __init__(self, X, Y, k_fold):  
        self.X = np.array(X)
        self.Y = np.array(Y)
        self.skf = StratifiedKFold(n_splits = k_fold)
        self._kfold = list(self.skf.split(self.X, self.Y))

    def _use_ith_fold(self, kth):
        if kth >= len(self._kfold):
            raise ValueError("The kth-fold is already exhausted.")
        kth_fold = self._kfold[kth]
        return kth_fold


          