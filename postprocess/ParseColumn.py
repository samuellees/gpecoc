#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2018-05-08 10:17
@author: Eric

# 
"""
import os

from Base import Parser
from utils.stack import Stack
import Configurations as Configs
from utils.dirtools import check_folder, del_dir_tree


class ParseColumn(object):

    __root__ = None
    __currentdataset__ = None
    __currentdir__ = None
    __currentinfile__ = None
    __currentoutfile__ = None
    __rdfile__ = None # file under reading
    __wrtfile__ = None # file under writing
    __datasets__ = None
    __dirs__ = None
    __files__ = None
    __genidset__ = None

    def __init__(self):
        pass

    def setRoot(self, r_path):
        self.__root__ = r_path

    def setGenids(self, genids):
        self.__genidset__ = set()
        for gid in genids:
            self.__genidset__.add(gid)

    # dataset
    def readDatasets(self):
        self.__datasets__ = Stack()
        flist = os.listdir(self.__root__)
        for fname in flist:
            fpath = os.path.join(self.__root__, fname)
            if (os.path.isdir(fpath)) and ('-v' in fname):
                self.__datasets__.push(fpath)

    def hasNextDataset(self):
        return not self.__datasets__.isEmpty()

    def nextDataset(self):
        if not self.__datasets__.isEmpty():
            self.__currentdataset__ = self.__datasets__.pop()

    # dir
    def readDirs(self):
        self.__dirs__ = Stack()
        flist = os.listdir(self.__currentdataset__)
        for fname in flist:
            fpath = os.path.join(self.__currentdataset__, fname)
            if (os.path.isdir(fpath)) and ('hamm' in fname):
                self.__dirs__.push(fpath)

    def hasNextDir(self):
        return not self.__dirs__.isEmpty()

    def nextDir(self):
        if not self.__dirs__.isEmpty():
            self.__currentdir__ = self.__dirs__.pop()

    # file
    def readFiles(self):
        Parser.files = Stack()
        flist = os.listdir(self.__currentdir__)
        for fname in flist:
            fpath = os.path.join(self.__currentdir__, fname)
            if 'Gen' in fname and int(fname.split('.')[1]) in self.__genidset__:
                self.__files__.push(fpath)

    def hasNextFile(self):
        return not self.__files__.isEmpty()

    def nextFile(self):
        if not self.__files__.isEmpty():
            self.__currentinfile__ = self.__files__.pop()

    def openFileReader(self):
        self.__rdfile__ = open(self.__currentinfile__, 'r')

    def openFileWriter(self):
        self.__wrtfile__ = open(self.__currentoutfile__, 'wb+')

    def writeLine(self, line):
        self.__wrtfile__.writelines(line+'\n')
        self.__wrtfile__.flush()

    def closeReader(self):
        self.__rdfile__.close()
        self.__rdfile__ = None

    def closeWriter(self):
        self.__wrtfile__.flush()
        self.__wrtfile__.close()
        self.__wrtfile__ = None

    def setInFile(self, genid):
        self.__currentinfile__ = os.path.join(self.__currentdir__, "Gen."+str(genid)+".gpecoc")

    def setOutFile(self):
        if '\\' in self.__currentdataset__:
            datasetName = self.__currentdataset__.split("\\")[-1].split('-')[0]
        else:
            datasetName = self.__currentdataset__.split("/")[-1].split('-')[0]
        fpath = os.path.join(self.__root__, 'a_s'+Configs.version)
        fpath = os.path.join(fpath, 'a_Column')
        check_folder(fpath)
        self.__currentoutfile__ = os.path.join(fpath, datasetName)
        del_dir_tree(self.__currentoutfile__)

    def parse_column(self):
        nextgenid = -1
        string = ''

        while True:
            nextgenid = nextgenid+1
            if nextgenid >= Configs.generations:
                break
            if (nextgenid+1) not in self.__genidset__:
                continue

            self.setInFile(nextgenid)
            self.openFileReader()
            reader = self.__rdfile__

            # origin column
            reader.readline()
            line = reader.readline()
            origin_column = (len(line)-2)/3

            # add columns
            add_columns =0
            line = reader.readline()
            while "2:" not in line:
                if "Add one column:" in line:
                    add_columns+=1
                line = reader.readline()

            column = origin_column + add_columns
            string = string + '%d\t' % column
            self.closeReader()
        return string


    def parse(self, root, genids):
        self.setGenids(genids)
        self.setRoot(root)
        self.readDatasets()
        while(self.hasNextDataset()):
            self.nextDataset()
            self.readDirs()
            self.setOutFile()
            self.openFileWriter()
            while(self.hasNextDir()):
                self.nextDir()
                line = self.parse_column()
                self.writeLine(line)
            self.closeWriter()

