# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 13:32:34 2016

@author: hanrui
"""
import copy

import numpy as np


import utils.LegalityCheckers as LC
from ecoc.OperationFuncs import *

nodeType = {"TERMINAL": 0, "NONTERMINAL": 1}

"""
# convert Tree to EcocMatrix, 
# and then use LegalityCheckers to delete illegal columns
"""
def getMatrixDirectly_and_feature(ind):
    import GTree
    arrays = []
    feature = []
    prefixnum = 0
    features = gol.get_val("feature_method")
    classes = gol.get_val("classes")
    MaxDeapth = gol.get_val("maxDeap")
    New_Ind = GTree.GTreeGP()
    for i in xrange(len(ind.nodes_list)):
        if ind.nodes_list[i].getType() == nodeType["NONTERMINAL"]:
            prefixnum = prefixnum + 1
            for j in xrange(0, len(classes)):
                locals()[classes[j]] = classes[j]

            New_Ind.setRoot(ind.nodes_list[i])
            array = eval(New_Ind.getCompiledCode())
            arrays.append(list(array))
        else:
            arrays.append(list(ind.nodes_list[i].getData()))
      
    for i in xrange(len(arrays)):
        for j in xrange(len(features)):
            if features[j] in arrays[i]:
                feature.append(features[j])
                arrays[i].remove(features[j])
                
    pre_order = ind.nodes_list
    cem = CEM(MaxDeapth, classes, prefixnum)
    cem.construct_tree(pre_order, arrays)
    ecocMatrix = np.array(cem.getMatrix())
    feature = np.array(feature)
    

    #1.There being a column with all 0 or 1 or -1
    deletes = LC.zeroColumn(ecocMatrix)
    ecocMatrix = np.delete(ecocMatrix, deletes, axis=1)
    feature = np.delete(feature, deletes, axis=0)

    #2.There being a column that is lack of 1 or -1
    deletes = LC.onlyOneColumn(ecocMatrix)
    ecocMatrix = np.delete(ecocMatrix, deletes, axis=1)
    feature = np.delete(feature, deletes, axis=0)
    
    #3.Two columns having the same numbers
    deletes = LC.sameColumns(ecocMatrix, feature)
    ecocMatrix = np.delete(ecocMatrix, deletes, axis=1)
    feature = np.delete(feature, deletes, axis=0)
    
    #4.Two columns having the opposite numbers
    deletes = LC.opstColumns(ecocMatrix, feature)
    ecocMatrix = np.delete(ecocMatrix, deletes, axis=1)
    feature = np.delete(feature, deletes, axis=0)
    #############################
    return ecocMatrix, feature


def getMatrixDirectly(ind):
    ind_info = getMatrixDirectly_and_feature(ind)
    return ind_info[0]


"""
# base class for converting tree to ecocmatrix
"""
class CEM:

    max = 1000000
     
    def __init__(self, d, c, p):  
        
        self.matrix = []
        self.cengl = []
        self.cengh = []
        self.temp = []
        self.prefixNumcount = 0
        self.di = -1
        self.features = []
        self.delete_flag = False
        
        # INIT MATRIX
        self.classes = c
        self.m = []
        for i in range(len(self.classes)):  
            self.m.append(0)
        self.matrix.append(self.m) 
        
        # INIT 
        self.deap = d
        for i in range(self.deap):
            self.temp.append(0)
            
        # cengl从第几个非终端节点开始； cengh到达终点时跨越了几步
        self.prefixNum = p
        for i in range(self.prefixNum): 
            self.cengl.append(self.max)   
            self.cengh.append(self.max)
            self.features.append(" ")
     
            
    def construct_tree(self, pre_order, brrays):
        
        if len(pre_order) == 0 :
            return None
            
        if pre_order[0].getType() == nodeType["NONTERMINAL"]:
            self.prefixNumcount = self.prefixNumcount + 1
            self.di = self.di + 1
            #print "this is not a terminal"
            #print brrays[0]
            
            # init current matrix
            for i in range(len(self.classes)):  
                self.m[i] = 0
            
            self.cengl[self.di] = self.di
            self.cengh[self.di] = self.cengl[self.di] + 1
            
            self.temp[self.di] = self.prefixNumcount
            
            # +1
            for s in brrays[1]:
                self.m[self.classes.index(s)] = 1
            # add list
            x = copy.deepcopy(self.m)
            self.matrix.append(x)
              
            '''
            print "---------------left---------------------------------------------"
            print self.matrix
            print "----------------------------------------------------------------"
            '''
            
            # add step for left
            for i in range(len(self.cengh)):
                self.cengh[i] = self.cengh[i] + 1
    
            self.construct_tree(pre_order[1:], brrays[1:])
            
            #print "left is over, right is begining"
            
            # add step for right
            for i in range(len(self.cengh)):
                self.cengh[i] = self.cengh[i] + 1
    
            '''
            print "---------------right-------------------------------------------"
            print "brrays[(cengh[di] - cengl[di] - 1)"
            print brrays[(self.cengh[self.di] - self.cengl[self.di] - 1)]
            '''
            # -1
            for s in brrays[(self.cengh[self.di] - self.cengl[self.di] - 1)]:
                if self.matrix[self.temp[self.di]][self.classes.index(s)] == 1:
                    self.matrix[self.temp[self.di]][self.classes.index(s)] = 0
                else:
                    self.matrix[self.temp[self.di]][self.classes.index(s)] = -1
  
            '''
            print self.matrix
            print "---------------------------------------------------------------"
            
            print self.cengh[self.di] - self.cengl[self.di]
            '''
            self.construct_tree(pre_order[(self.cengh[self.di] - self.cengl[self.di] - 1):], brrays[(self.cengh[self.di] - self.cengl[self.di] - 1):])
            
            # reback to the last prefix dot
            self.di = self.di - 1
    
        #else:
            
            #print "this is a terminal"
            #print brrays[0]
        return ;

        
    def display(self):
        del self.matrix[0]
        for each_list in self.matrix:  
            print (each_list)  
      
            
    def getMatrix(self):
        if self.delete_flag == False:
            del self.matrix[0]
            self.delete_flag = True
        return map(list, zip(*self.matrix))
          
