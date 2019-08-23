# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 11:40:22 2017

@author: Shone
"""

import os


# if the dir not exist, then create it 
def check_folder(_path):
    if not os.path.isdir(_path):
        os.makedirs(_path)


# delete the dir tree
def del_dir_tree(path):
    if os.path.isfile(path):
        try:
            os.remove(path)
        except Exception, e:
            #pass
            print e
    elif os.path.isdir(path):
        for item in os.listdir(path):
            itempath = os.path.join(path, item)
            del_dir_tree(itempath)
        try:
            os.rmdir(path)   # 
        except Exception, e:
            #pass
            print e