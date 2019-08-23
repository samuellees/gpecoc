#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from sklearn.metrics.pairwise import euclidean_distances

import Configurations as Configs
import utils.gol as gol
import ecoc.CallbackFuncs as CB
import ecoc.EvaluateMethods as EM
import ecoc.Initializator as Initializator
import gp.TreeMatrixConvertor as TMConvertor
from ecoc.ConnectClassifier import ConnectClassifier
from ecoc.OperationFuncs import *
from gp import Consts
from gp import GSimpleGA, GTree


def main_run():
    ##########################################
    # variables preparation
    ##########################################
    Initializator.init_gol()
    gol.set_val("aimFolder", Configs.aimFolder)
    gol.set_val("dataName", Configs.dataName)
    Initializator.init_all()
    classes = gol.get_val("classes")
    maxDeap = gol.get_val("maxDeap")
    growMethod = gol.get_val("growMethod")
    generations = gol.get_val("generations")
    crossoverRate = gol.get_val("crossoverRate")
    mutationRate = gol.get_val("mutationRate")
    populationSize = gol.get_val("populationSize")
    freq_Stats = gol.get_val("freq_stats")
    Train_X = gol.get_val("Train_X")
    Train_Y = gol.get_val("Train_Y")
    validation_X = gol.get_val("validation_X")
    validation_Y = gol.get_val("validation_Y")
    Test_X = gol.get_val("Test_X")
    Test_Y = gol.get_val("Test_Y")
    sel_features = gol.get_val("sel_features")
    ##########################################

    genome = GTree.GTreeGP()
    genome.setParams(max_depth=maxDeap, method=growMethod)
    genome.evaluator += EM.eval_func_fscore


    ga = GSimpleGA.GSimpleGA(genome)
    ga.setParams(gp_terminals=classes, gp_function_prefix="Operation")
    ga.setMinimax(Consts.minimaxType["maximize"])
    ga.setGenerations(generations)
    ga.setCrossoverRate(crossoverRate)
    ga.setMutationRate(mutationRate)
    ga.setPopulationSize(populationSize)
    ga.setElitismReplacement(1)
    #ga.stepCallback.set(CB.printIndividuals_callback)
    ga.stepCallback += CB.checkAncients_callback
    ga.stepCallback += CB.logResultEveryGen_callback
    ga.stepCallback += CB.delogPopulation_callback
    ga.stepCallback += CB.logMiddleInfo_callback
    ga.stepCallback += CB.debug_callback
    
    print "------------------------------------------------------"
    
    ga(freq_stats=freq_Stats)
    best = ga.bestIndividual()

    #change the display_flag to display test labels and predict labels
    FinalMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(best)
    cc = ConnectClassifier(features_used_list, sel_features , FinalMatrix)
    finalScore, finalAccuracy, infos_evaluations = cc.FinalTrainAndTest(Train_X, Train_Y, validation_X, validation_Y, Test_X, Test_Y)

    # euddist
    num_class = len(classes)
    num_cols = FinalMatrix.shape[1]
    _dist = euclidean_distances(FinalMatrix, FinalMatrix)/np.sqrt(num_cols)
    dist = np.sum(_dist)/2/(num_class*(num_class-1))

    infos_evaluations.insert(len(infos_evaluations), "---------test------------")
    infos_evaluations.insert(len(infos_evaluations), "fscore: %f" % finalScore)
    infos_evaluations.insert(len(infos_evaluations), "accuracy: %f" % finalAccuracy)
    infos_evaluations.insert(len(infos_evaluations), "dist: %f" % dist)

    for text in infos_evaluations:
        print text


if __name__ == "__main__":
    main_run()
