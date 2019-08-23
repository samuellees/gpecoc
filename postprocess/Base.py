# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 17:41:14 2017

@author: root
"""

import os

from abc import ABCMeta, abstractmethod
from utils.stack import Stack


class Parser(object):

    __metaclass__ = ABCMeta
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

    @abstractmethod
    def setInFile(self): pass

    @abstractmethod
    def parseFile(self): pass

    @abstractmethod
    def setOutFile(self): pass

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
				self.setInFile()
				self.openFileReader()
				line = self.parseFile()
				self.writeLine(line)
				self.closeReader()
			self.closeWriter()



