#!/usr/bin/env python
# -*- coding: utf-8 -*-
import getpass
'''
# Basic Configuration

version: 
  Must be changed every time you run the program, 
  or data will be over written and become useless.

dataName: 
  Used for debuging.

aimFolder: 
  Used for debuging.

generations: 
  #generations for Genetic Programming.

populationSize: 
  populationSize for Genetic Programming.

n_jobs: 
  setted to 1, don't change that.

freq_stats: 
  setted to 1, don't change that.

n_neighbors: 
  hyper params for KNN, it will be used when you choose KNN as base classifier.

mutationRate, crossoverRate:
  mutationRate, crossoverRate for Genetic Programming.

growMethod:
  setted to "ramped" as mentioned in the paper, don't change that.

root_path:
  The root path for your project, the suggested path is D:\gpecoc.
  NOTE: permission ISSUES will occur is you put the project in 'C' disk(system disk).

'''


version = "8.80"

dataName = "vertebral"

aimFolder = "hamm1"

generations = 60

populationSize = 10

n_jobs = 1

freq_stats = 1

n_neighbors = 3

crossoverRate = 1

mutationRate = 0.1

growMethod = "ramped"

root_path = "D:\parallel\gpecoc1"
