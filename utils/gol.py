# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:24:35 2017

To save some golbal variables, like training set and so on.

@author: Shone
"""



def _init():#初始化
    global _global_dict
    _global_dict = {}


def set_val(key,value):
    _global_dict[key] = value


def get_val(key,defValue=None):
    try:
        return _global_dict[key]
    except KeyError:
        return defValue