# -*- coding: utf-8 -*-

import os
import gol
import sys
import numpy as np

from utils.dirtools import check_folder
from gp.TreeMatrixConvertor import getMatrixDirectly_and_feature


def deprint( tag , value ):
    print "===%s====start====" %tag
    print value
    print "===%s=====end=====" %tag
    

def deprint_string(string):
    #None
    print string
    
def decache(string):
    #None
    sys.stdout.write(string)
    sys.stdout.flush()
    
    
def logPopulations( genid, pop):
    filedir = gol.get_val("root_path")
    filedir = os.path.join(filedir, 'Results/' + gol.get_val("version"))
    filedir = os.path.join(filedir, gol.get_val("dataName")+"-v"+gol.get_val("version"))
    filedir = os.path.join(filedir, gol.get_val("aimFolder"))
    check_folder(filedir)
    filedir = os.path.join(filedir, "Gen."+str(genid)+".gpecoc")

    f = file(filedir, 'w+')
    i = 1
    for ind in pop:
        Matrix , feature_list = getMatrixDirectly_and_feature(ind)
        feature_method_index = gol.get_val("feature_method_index")
        feature_index_list = list( feature_method_index[method] for method in feature_list)
        Matrix = np.ndarray.tolist(Matrix)
        Matrix.insert(0,feature_index_list)
        Matrix = np.array(Matrix)
        
        f.write(str(i))
        f.write(":")
        f.write('\n')
        f.write(str(Matrix))
        f.write('\n')
        f.write("f-score:")
        f.write(str(ind.fscore))
        f.write('\n')
        f.write("accuracy:")
        f.write(str(ind.accuracy))
        f.write('\n')
        f.write('--------------------------------------------------------------------------------------------\n')
        
        for text in ind.infos_evaluation:
            f.write(str(text))
            f.write('\n')
        f.write('\n\n\n\n\n')
        i=i+1
    f.close()
    
    
def logMiddle(genid, fitness, filename):
    filedir = gol.get_val("root_path")
    filedir = os.path.join(filedir, 'Results/' + gol.get_val("version"))
    filedir = os.path.join(filedir, gol.get_val("dataName")+"-v"+gol.get_val("version"))
    filedir = os.path.join(filedir, gol.get_val("aimFolder"))
    check_folder(filedir)
    filedir = os.path.join(filedir, filename)
    
    f = file(filedir, 'a')
    i = genid+1
    f.write(str(i))
    f.write(":")
    f.write('\n')
    f.write("fitness:")
    f.write(str(fitness))
    f.write('\n')
    f.write('\n')
    f.close()