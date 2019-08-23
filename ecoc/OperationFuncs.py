#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
# there are some operation methods
'''

from utils import gol

# there is hard code in Crossover module
# self-defined operation method, prefix "Operation" is necessary
def Operation_F1(a, b): 
    features = gol.get_val("feature_method")
    a = list(a)
    b = list(b)
    for i in xrange(len(features)):
        if features[i] in a:
            a.remove(features[i])
        if features[i] in b:
            b.remove(features[i])
    sets = set(a)|set(b)
    function = "fclassify"
    result = list(sets)  
    result.insert(0, function)
    return result
    
def Operation_F2(a, b): 
    features = gol.get_val("feature_method")
    a = list(a)
    b = list(b)
    for i in xrange(len(features)):
        if features[i] in a:
            a.remove(features[i])
        if features[i] in b:
            b.remove(features[i])
    
    sets = set(a)|set(b)
    function = "svmrfe"
    result = list(sets)  
    result.insert(0, function)
    return result
    
def Operation_F3(a, b): 
    features = gol.get_val("feature_method")
    a = list(a)
    b = list(b)
    for i in xrange(len(features)):
        if features[i] in a:
            a.remove(features[i])
        if features[i] in b:
            b.remove(features[i])
    
    sets = set(a)|set(b)
    function = "forest"
    result = list(sets)  
    result.insert(0, function)
    return result  
    
def Operation_F4(a, b): 
    features = gol.get_val("feature_method")
    a = list(a)
    b = list(b)
    for i in xrange(len(features)):
        if features[i] in a:
            a.remove(features[i])
        if features[i] in b:
            b.remove(features[i])
    
    sets = set(a)|set(b)
    function = "bsswss"
    result = list(sets)  
    result.insert(0, function)
    return result  