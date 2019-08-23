# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 19:41:12 2017

@author: Shone
"""


import os

from Base import Parser
import Configurations as Configs
from utils.dirtools import check_folder, del_dir_tree


class ParserAHD(Parser):
    
    def __init__(self): None
    
    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "Result")
    
    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_AHD')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)
        
    def parseFile(self):
        nextgenid = -1
        reader = self.__rdfile__
        string = ''
        
        while True:
            line = reader.readline()
            if 'Gen. ' not in line: 
                continue
            
            nextgenid = nextgenid+1
            if nextgenid >= Configs.generations:
                break
            
            if (nextgenid+1) not in self.__genidset__:
                continue
            
            ss = line.split()
            ahd = float(ss[6])
            string = string + '%.4f\t' % ahd
            
        return string


class ParserFitness(Parser):
    
    def __init__(self): None
    
    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "Result")
    
    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_Fitness')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)
        
    def parseFile(self):
        nextgenid = -1
        reader = self.__rdfile__
        string = ''
        
        while True:
            line = reader.readline()
            if 'Gen. ' not in line: 
                continue
            
            nextgenid = nextgenid+1
            if nextgenid >= Configs.generations:
                break
            
            if (nextgenid+1) not in self.__genidset__:
                continue
            
            ss = line.split()
            fitness = float(ss[2])
            string = string + '%.4f\t' % fitness
            
        return string


class ParserTestAcc(Parser):
    
    def __init__(self): None
    
    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "AAAAA")
    
    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_TestAcc')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)
        
    def parseFile(self):
        nextgenid = 0
        reader = self.__rdfile__
        string = ''
        
        while True:
            reader.readline()
            line = reader.readline()
            reader.readline()
            
            nextgenid = nextgenid+1
            if nextgenid > Configs.generations:
                break
            
            if nextgenid not in self.__genidset__:
                continue
            
            ss = line.split(':')
            accuracy = float(ss[1])
            string = string + '%.4f\t' % accuracy
            
        return string


class ParserTestFscore(Parser):

    def __init__(self): None

    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "AAAAAfscore")

    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_TestFscore')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)

    def parseFile(self):
        nextgenid = 0
        reader = self.__rdfile__
        string = ''

        while True:
            reader.readline()
            line = reader.readline()
            reader.readline()

            nextgenid = nextgenid+1
            if nextgenid > Configs.generations:
                break

            if nextgenid not in self.__genidset__:
                continue

            ss = line.split(':')
            fscore = float(ss[1])
            string = string + '%.4f\t' % fscore

        return string


class ParserTrainAcc(Parser):
    
    def __init__(self): None
    
    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "Result")
    
    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_TrainAcc')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)
        
    def parseFile(self):
        nextgenid = -1
        reader = self.__rdfile__
        string = ''
        
        while True:
            line = reader.readline()
            if 'Gen. ' not in line: 
                continue
            
            nextgenid = nextgenid+1
            if nextgenid >= Configs.generations:
                break
            
            if (nextgenid+1) not in self.__genidset__:
                continue
            
            ss = line.split()
            acc = float(ss[7])
            string = string + '%.4f\t' % acc
            
        return string


class ParserTrainFscore(Parser):
    
    def __init__(self):
        pass
    
    def setInFile(self):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "Result")
    
    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_TrainFscore')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)
        
    def parseFile(self):
        nextgenid = -1
        reader = self.__rdfile__
        string = ''
        
        while True:
            line = reader.readline()
            if 'Gen. ' not in line: 
                continue
            
            nextgenid = nextgenid+1
            if nextgenid >= Configs.generations:
                break
            
            if (nextgenid+1) not in self.__genidset__:
                continue
            
            ss = line.split()
            fscore = float(ss[5])
            string = string + '%.4f\t' % fscore
            
        return string

