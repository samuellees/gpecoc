# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 11:16:17 2017

# Data loading modules

@author: Shone
"""

import numpy as np
from utils import delog
from numpy import double, loadtxt
from sklearn import preprocessing


'''
# load lables dictionary in form of 'a,b,c,d' , from files directly.
'''
def loadClasses(trainFile, validationFile, testFile):
    trainStrY = loadStrLabel(trainFile)
    validationStrY = loadStrLabel(validationFile)
    testStrY = loadStrLabel(testFile)
    _uniqueStrY =  np.array([e for e in np.unique(trainStrY) if e in np.unique(testStrY)])
    uniqueStrY = np.array([e for e in np.unique(_uniqueStrY) if e in np.unique(validationStrY)])
    np.sort(uniqueStrY)
    classes = list()
    for i in xrange(len(uniqueStrY)):
        ascii = ord('A') + i
        if(ascii>90):
            ascii += 6
        classes.append(chr(ascii))
    return classes
    

'''
# load scaled data, python list, a row is a simple
'''
def loadDataset(trainFile, validationFile, testFile):
    trainX, validationX, testX = loadScaledFeatures(trainFile, validationFile, testFile)
    trainStrY = loadStrLabel(trainFile)
    validationStrY = loadStrLabel(validationFile)
    testStrY = loadStrLabel(testFile)
    length = len(trainStrY) + len(validationStrY) + len(testStrY)
    class_num = max(len(np.unique(trainStrY)),len(np.unique(validationStrY)),len(np.unique(testStrY)))
    trainX, trainStrY, testX, testStrY = labelChecker(trainX, trainStrY, testX, testStrY)
    trainX, trainStrY, validationX, validationStrY = labelChecker(trainX, trainStrY, validationX, validationStrY)
    validationX, validationStrY, testX, testStrY = labelChecker(validationX, validationStrY, testX, testStrY)
    trainY, validationY, testY = str2alpha(trainStrY, validationStrY, testStrY)

    return trainX, trainY, validationX, validationY, testX, testY

'''
# load feature data, python list, z-score, a row is a simple
'''
def loadScaledFeatures(trainFile, validationFile, testFile):
    raw_trainX= np.transpose(np.loadtxt(trainFile,skiprows=1,dtype=double,ndmin=2,delimiter=','))
    raw_validationX= np.transpose(np.loadtxt(validationFile,skiprows=1,dtype=double,ndmin=2,delimiter=','))
    raw_testX= np.transpose(np.loadtxt(testFile,skiprows=1,dtype=double,ndmin=2,delimiter=','))
    scaler = preprocessing.StandardScaler().fit(raw_trainX)
    trainX = scaler.transform(raw_trainX)
    validationX = scaler.transform(raw_validationX)
    testX = scaler.transform(raw_testX)
    return trainX, validationX, testX



'''
# read lable, return list of string
'''
def loadStrLabel(filename):
        file_ = open(filename, "r")
        line = file_.readline()
        line = line.strip() #delete the last character
        labels = line.split(',')
        file_.close()
        return labels


'''
# check and do some ajustment on the lables in training set and testing set.
'''
def labelChecker(trainX, trainStrY, testX, testStrY):
    trainX = np.array(trainX)
    trainStrY = np.array(trainStrY)
    testX = np.array(testX)
    testStrY = np.array(testStrY)
    
    train_label_unique = np.unique(trainStrY)
    test_label_unique = np.unique(testStrY)
    useful_label_unique = np.array([e for e in train_label_unique if e in test_label_unique ])
    if len(train_label_unique) > len(useful_label_unique) :
        useless_label = np.array([e for e in  (set(train_label_unique) - set(useful_label_unique))])
        useless_index = np.array([i for i in xrange(len(trainStrY)) 
                                                if trainStrY[i] in useless_label])                                     
        #delete rows
        trainStrY = np.delete(trainStrY,useless_index,0)
        trainX = np.delete(trainX,useless_index,0)
        delog.deprint_string("Useless labels in training samples:")
        delog.deprint_string( useless_label)
        delog.deprint_string("Data about labels above have been removed from training data.")
    if len(test_label_unique) > len(useful_label_unique) :
        useless_label = np.array([e for e in (set(test_label_unique) - set(useful_label_unique))])
        useless_index = np.array([i for i in xrange(len(testStrY)) 
                                                if testStrY[i] in useless_label])                                    
        #delete rows
        testStrY = np.delete(testStrY,useless_index,0)
        testX = np.delete(testX,useless_index,0)
        delog.deprint_string("Useless labels in testing samples:")
        delog.deprint_string(useless_label)
        delog.deprint_string("Data about labels above have been removed from testing data.")
    if len(train_label_unique) == len(useful_label_unique) and len(test_label_unique) == len(useful_label_unique):
        delog.deprint_string("No error in training data and testing data.")
    return trainX, trainStrY, testX, testStrY
        

def str2alpha(trainStrY, validationStrY, testStrY):
    uniqueStrY = np.unique(trainStrY)
    np.sort(uniqueStrY)
    str2alpha_dict = dict()
    for i in xrange(len(uniqueStrY)):
        ascii = ord('A') + i
        if(ascii>90):
            ascii += 6
        str2alpha_dict[uniqueStrY[i]] = chr(ascii)
    
    trainY = list()
    validationY = list()
    testY = list()
    for i in xrange(len(trainStrY)):
        trainY.append(str2alpha_dict[trainStrY[i]])
    for k in xrange(len(validationStrY)):
        validationY.append(str2alpha_dict[validationStrY[k]])
    for j in xrange(len(testStrY)):
        testY.append(str2alpha_dict[testStrY[j]])
    
    trainY, validationY, testY = np.array(trainY),np.array(validationY),np.array(testY)
    return trainY, validationY, testY
