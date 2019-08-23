#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys

sys.path.append(os.getcwd)

"""
:mod:`pyevolve` -- the main pyevolve namespace
================================================================

This is the main module of the pyevolve, every other module
is above this namespace, for example, to import :mod:`Mutators`:

"""


__version__ =  '0.6rc1'
__author__ =  'Christian S. Perone'

import sys
import logging
import gp.Consts as Consts

if sys.version_info[:2] < (2, 5):
   raise Exception("Python 2.5+ required, the version %s was found on your system !" % (sys.version_info[:2],))

del sys

def logEnable(filename=Consts.CDefLogFile, level=Consts.CDefLevel):
   """ Enable the log system for pyevolve

   :param filename: the log filename
   :param level: the debugging level

   Example:
      >>> pyevolve.logEnable()

   """
   logging.basicConfig(level=level,
                    format='%(asctime)s [%(module)s:%(funcName)s:%(lineno)d] %(levelname)s %(message)s',
                    filename=filename,
                    filemode='w')
   logging.info("Pyevolve v.%s, the log was enabled by user.", __version__)