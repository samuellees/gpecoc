# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:57:31 2017

@author: Shone
"""

import random

import copy
import numpy as np

import utils.LegalityCheckers as LCheckers
import gp.TreeMatrixConvertor as TMConvertor
from ConnectClassifier import ConnectClassifier as CC
from gp import Consts
from gp import Util, GTree
from gp import GTreeNode
from utils import gol


def debug_callback(gp_engine):
    genes = gol.get_val("genebank").genes
    genid = gp_engine.getCurrentGeneration()
    None


# log the best accuracy of every generation into file AAAAA
# log the best fscore in every generation into file AAAAAfscore
# log the best accuracy without local improvement into file BestAcc_no_impro
# log the best fscore without local improvement into file BestFscore_no_impro
def logMiddleInfo_callback(gp_engine):
    Train_X = gol.get_val("Train_X")
    Train_Y = gol.get_val("Train_Y")
    validation_X = gol.get_val("validation_X")
    validation_Y = gol.get_val("validation_Y")
    Test_X = gol.get_val("Test_X")
    Test_Y = gol.get_val("Test_Y")
    sel_features = gol.get_val("sel_features")

    import sys
    from utils import delog
    sys.stdout.write("logMiddleInfo...")
    genid = gp_engine.getCurrentGeneration()

    best = gp_engine.bestIndividual()
    FinalMatrix, features_used_list = TMConvertor.getMatrixDirectly_and_feature(best)

    # result with local improvemtent
    cc = CC(features_used_list, sel_features, FinalMatrix)
    finalScore, finalAccuracy, infos_evaluations = cc.FinalTrainAndTest(Train_X, Train_Y, validation_X, validation_Y, Test_X, Test_Y)
    delog.logMiddle(genid, finalAccuracy, "AAAAA")
    delog.logMiddle(genid, finalScore, "AAAAAfscore")

    #  result without local improvemtent
    cc = CC(features_used_list, sel_features, FinalMatrix)
    cc.TrainAndTest_withoutlocalimp(Train_X, Train_Y, validation_X, validation_Y)
    _finalScore, _finalAccuracy = cc.FinalTrainAndTest_withoutlocalimp(Train_X, Train_Y, validation_X, validation_Y, Test_X, Test_Y)
    delog.logMiddle(genid, _finalAccuracy, "BestAcc_no_impro")
    delog.logMiddle(genid, _finalScore, "BestFscore_no_impro")

    sys.stdout.write("over\n")
    sys.stdout.flush()


def delogPopulation_callback(gp_engine):
    from utils import delog
    pop = gp_engine.getPopulation()
    genid = gp_engine.getCurrentGeneration()
    delog.logPopulations(genid,pop)


def logResultEveryGen_callback(gp_engine):
    if gp_engine.getCurrentGeneration() ==0:
        print "="*65
        format_str = 'Gen' + ' '*12 + '%%-8s  %%-8s  %%-8%s %%-10%s   %%-10%s   %%-10%s'
        print( (format_str % ('s', 's', 's', 's')) % ('Max', 'Min', 'Avg', 'Best-Fscore', 'Best-Hamdist', 'Best-Accuracy'))
    np.set_printoptions(threshold='nan') 
    # do in every generation
    best = gp_engine.getPopulation().bestRaw()
    bestMatrix , feature_list = TMConvertor.getMatrixDirectly_and_feature(best)
    feature_method_index = gol.get_val("feature_method_index")
    feature_index_list = list(feature_method_index[method] for method in feature_list)
    bestMatrix = np.ndarray.tolist(bestMatrix)
    bestMatrix.insert(0,feature_index_list)
    print np.array(bestMatrix)




def checkAncients_callback(gp_engine):
    if gp_engine.getCurrentGeneration() != 0: return
    from utils import delog
    delog.decache("check first Gen...")

    begin = 0
    end = gol.get_val("populationSize")
    classes = gol.get_val("classes")
    population = gp_engine.getPopulation()
    for i in xrange(begin, end):
        genome = population[i]
        max_depth = genome.getParam("max_depth", None)

        #illegal?
        ecocMatrix, feature_list = TMConvertor.getMatrixDirectly_and_feature(genome)
        Illegal = False
        if LCheckers.tooLittleColumn(ecocMatrix):
            Illegal = True
        elif LCheckers.tooMuchColumn(ecocMatrix):
            Illegal = True
        # 2. if any class not included in the terminal nodes.
        else:
            labels = set(classes)
            for i in genome.nodes_list:
                if i.isLeaf():
                    labels = labels - set(i.getData())
            labels = list(labels)
            if len(labels) > 0:
                Illegal = True

        if max_depth is None:
            Util.raiseException("You must specify the max_depth genome parameter !", ValueError)
        if max_depth < 0:
            Util.raiseException("The max_depth must be >= 1, if you want to use GTreeGPMutatorSubtree crossover !", ValueError)

        while Illegal==True:
            new_genome = copy.deepcopy(genome)
            node = new_genome.getRandomNode()
            assert node is not None
            depth = new_genome.getNodeDepth(node)
            node_parent = node.getParent()
            root_subtree = GTreeNode.buildGTreeGPGrow(gp_engine, 0, max_depth - depth)
            if node_parent is None:
                new_genome.setRoot(root_subtree)
            else:
                root_subtree.setParent(node_parent)
                node_parent.replaceChild(node, root_subtree)
            new_genome.processNodes()

            # illegal ? 
            # Actually, case #1 and case #2 may not happen
            Illegal = False
            ecocMatrix, feature_list = TMConvertor.getMatrixDirectly_and_feature(new_genome)

            # 1.The number of column is too little
            if LCheckers.tooLittleColumn(ecocMatrix):
                Illegal = True
            elif LCheckers.tooMuchColumn(ecocMatrix):
                Illegal = True
            # 2. if any class not included in the terminal nodes.
            else:
                labels = set(classes)
                for i in new_genome.nodes_list:
                    if i.isLeaf():
                        labels = labels - set(i.getData())
                labels = list(labels)
                if len(labels) > 0:
                    Illegal = True

            # apply the mutations
            if Illegal == False:
                genome.setRoot(new_genome.getRoot())
                genome.processNodes()

    #Update the scores of population
    delog.deprint_string( "over.")
    population.evaluate()
    population.sort()







numnum = 0
# self-defined printer
def printIndividuals_callback(gp_engine):
    import pydot
    global numnum
    New_Ind = GTree.GTreeGP()
    classes = gol.get_val("classes")
    numnum = numnum + 1
    begin = 0
    end = 20
    if gp_engine.getCurrentGeneration() != -1:
         population = gp_engine.getPopulation()
         graph = pydot.Dot(graph_type = "digraph")
         n = 0
         filename = 'Tree' + str(numnum) +'.jpg'
         for i in xrange(begin, end) :

             arrays = []
             ind = population[i]
             subg = pydot.Cluster("cluster_%d" % i, label="\"Ind. #%d - Score Raw/Fit.: %.4f/%.4f\"" % (i, ind.getRawScore(), ind.getFitnessScore()))
             count = n
             node_stack = []
             nodes_dict = {}
             tmp = None
             import __main__ as main_module

             for i in xrange(len(ind.nodes_list)):
                newnode = pydot.Node(str(count), style="filled")
                count += 1

                # color
                if ind.nodes_list[i].getType() == Consts.nodeType["TERMINAL"]:
                    newnode.set_color("lightblue2")
                else:
                   newnode.set_color("goldenrod2")

                # content of node
                if ind.nodes_list[i].getType() == Consts.nodeType["NONTERMINAL"]:
                    func = getattr(main_module, ind.nodes_list[i].getData())

                    if hasattr(func, "shape"):
                        newnode.set_shape(func.shape)

                    if hasattr(func, "representation"):
                        newnode.set_label(func.representation)
                    else:
                        for j in xrange(0, len(classes)):
                            locals()[classes[j]] = classes[j]

                        New_Ind.setRoot(ind.nodes_list[i])
                        array = eval(New_Ind.getCompiledCode())
                        newnode.set_label(str(array))
                    #if hasattr(func, "color"): newnode.set_color(func.color)
                else:
                    newnode.set_label(ind.nodes_list[i].getData())

                nodes_dict.update({ind.nodes_list[i]: newnode})
                graph.add_node(newnode)

             node_stack.append(ind.getRoot())
             while len(node_stack) > 0:
                tmp = node_stack.pop()

                parent = tmp.getParent()
                if parent is not None:
                   parent_node = nodes_dict[parent]
                   child_node  = nodes_dict[tmp]

                   newedge = pydot.Edge(parent_node, child_node)
                   graph.add_edge(newedge)

                rev_childs = tmp.getChilds()[:]
                rev_childs.reverse()
                node_stack.extend(rev_childs)
             n = count
             graph.add_subgraph(subg)
         graph.write(filename, prog='dot', format="jpeg")

