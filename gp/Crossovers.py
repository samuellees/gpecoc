"""

:mod:`Crossovers` -- crossover methods module
=====================================================================

In this module we have the genetic operators of crossover (or recombination) for each chromosome representation.

"""

import math
from random import randint as rand_randint, choice as rand_choice
from random import random as rand_random

import Util

nodeType = {"TERMINAL" : 0, "NONTERMINAL": 1}

#############################################################################
#################  GTreeGP Crossovers  ######################################
#############################################################################

def GTreeGPCrossoverSinglePoint(genome, **args):
   """ The crossover of the GTreeGP, Single Point for Genetic Programming

   ..note:: This crossover method creates offspring with restriction of the
            *max_depth* parameter.
   
   Accepts the *max_attempt* parameter, *max_depth* (required).   
   """
   #print "CrossoverAAAAAAAAAAA"
   sister = None
   brother = None

   gMom = args["mom"].clone()
   gDad = args["dad"].clone()

   gMom.resetStats()
   gDad.resetStats()

   max_depth   = gMom.getParam("max_depth", None)
   max_attempt = gMom.getParam("max_attempt", 15)

   if max_depth is None:
      Util.raiseException("You must specify the max_depth genome parameter !", ValueError)
      
   if max_depth < 0:
      Util.raiseException("The max_depth must be >= 1, if you want to use GTreeCrossoverSinglePointStrict crossover !", ValueError)

   momRandom = None
   dadRandom = None
   
   for i in xrange(max_attempt):

      dadRandom = gDad.getRandomNode()

      if   dadRandom.getType() == nodeType["TERMINAL"]:
         momRandom = gMom.getRandomNode(1)
      elif dadRandom.getType() == nodeType["NONTERMINAL"]:
         momRandom = gMom.getRandomNode(2)

      mD = gMom.getNodeDepth(momRandom)
      dD = gDad.getNodeDepth(dadRandom)

      # Two nodes are root
      if mD==0 and dD==0: continue

      mH = gMom.getNodeHeight(momRandom)
      if dD+mH > max_depth: continue

      dH = gDad.getNodeHeight(dadRandom)
      if mD+dH > max_depth: continue

      break

   if i==(max_attempt-1):
      assert gMom.getHeight() <= max_depth
      return (gMom, gDad)
   else:
      nodeMom, nodeDad = momRandom, dadRandom

   nodeMom_parent = nodeMom.getParent()
   nodeDad_parent = nodeDad.getParent()

   # Sister
   if args["count"] >= 1:
      sister = gMom
      nodeDad.setParent(nodeMom_parent)

      if nodeMom_parent is None:
         sister.setRoot(nodeDad)
      else:
         nodeMom_parent.replaceChild(nodeMom, nodeDad)
      sister.processNodes()
      assert sister.getHeight() <= max_depth

   # Brother
   if args["count"] == 2:
      brother = gDad
      nodeMom.setParent(nodeDad_parent)

      if nodeDad_parent is None:
         brother.setRoot(nodeMom)
      else:
         nodeDad_parent.replaceChild(nodeDad, nodeMom)
      brother.processNodes()
      assert brother.getHeight() <= max_depth

   return (sister, brother)



############
#self define
############

def _crossover_supplement_ecoc(raw_ind, coming_node):
   import numpy as np
   from utils import gol
   classes = gol.get_val("classes")
   leafs = raw_ind.getLeafs()
   datas = [_node_.getData() for _node_ in leafs]
   if len(classes)==len(np.unique(datas)):   return

   # when lack class
   import copy
   node_ind = copy.deepcopy(raw_ind)
   node_ind.setRoot(coming_node)
   node_ind.processNodes()
   leafs_coming = node_ind.getLeafs()

   leafs_old = copy.deepcopy(leafs)
   datas_old = [_node_.getData() for _node_ in leafs_old]
   for _node_ in leafs_coming: datas_old.remove(_node_.getData())

   class_lack = [cls for cls in classes if cls not in datas]

   # find nodes_redundancy, they will be replaced with lacks
   datas_temp_mid = []
   nodes_redundancy = []
   for _node_ in leafs_coming:
      if _node_.getData() in datas_old:   nodes_redundancy.append(_node_)
      elif _node_.getData() in datas_temp_mid:  nodes_redundancy.append(_node_)
      else: datas_temp_mid.append(_node_.getData())

   if len(nodes_redundancy)>0:
      for i in xrange(len(nodes_redundancy)):
         nodes_redundancy[i].setData(class_lack[0])
         class_lack.remove(class_lack[0])
         if len(class_lack)<=0: break

   # after replacing, it may still lack class
   if len(class_lack) > 0:
      import random
      from gp.GTreeNode import GTreeNodeGP
      max_depth = gol.get_val("maxDeap")
      nodes_grow_candidate = leafs_coming
      nodes_grow_candidate_new = nodes_grow_candidate

      # grow the tree to fill lacks
      while len(class_lack)>0:
         nodes_grow_candidate = nodes_grow_candidate_new
         nodes_grow_candidate_new = []
         for i in xrange(len(nodes_grow_candidate)):
            if raw_ind.getNodeDepth(nodes_grow_candidate[i]) < max_depth:
               _newnode_old = GTreeNodeGP(nodes_grow_candidate[i].getData(), node_type=nodeType['TERMINAL'], parent=nodes_grow_candidate[i])
               _newnode_lack = GTreeNodeGP(class_lack[0], node_type=nodeType['TERMINAL'], parent=nodes_grow_candidate[i])
               nodes_grow_candidate[i].addChild(_newnode_old)
               nodes_grow_candidate[i].addChild(_newnode_lack)
               nodes_grow_candidate[i].setType(nodeType['NONTERMINAL'])
               nodes_grow_candidate[i].setData('Operation_F' + str(random.randint(1,4)))
               class_lack.remove(class_lack[0])
               if len(class_lack)<=0: break
               nodes_grow_candidate_new.append(_newnode_old)
               nodes_grow_candidate_new.append(_newnode_lack)
         raw_ind.processNodes()
   raw_ind.processNodes()


def NewGTreeGPCrossover_ECOC(genome, **args):
   sister = None
   brother = None

   gMom = args["mom"].clone()
   gDad = args["dad"].clone()

   gMom.resetStats()
   gDad.resetStats()

   max_depth = gMom.getParam("max_depth", None)
   max_attempt = gMom.getParam("max_attempt", 15)

   if max_depth is None:
      Util.raiseException("You must specify the max_depth genome parameter !", ValueError)

   if max_depth < 0:
      Util.raiseException("The max_depth must be >= 1, if you want to use GTreeCrossoverSinglePointStrict crossover !",
                          ValueError)

   momRandom = None
   dadRandom = None

   for i in xrange(max_attempt):

      dadRandom = gDad.getRandomNode()

      if dadRandom.getType() == nodeType["TERMINAL"]:
         momRandom = gMom.getRandomNode(1)
      elif dadRandom.getType() == nodeType["NONTERMINAL"]:
         momRandom = gMom.getRandomNode(2)

      mD = gMom.getNodeDepth(momRandom)
      dD = gDad.getNodeDepth(dadRandom)

      # Two nodes are root
      if mD == 0 and dD == 0: continue

      mH = gMom.getNodeHeight(momRandom)
      if dD + mH > max_depth: continue

      dH = gDad.getNodeHeight(dadRandom)
      if mD + dH > max_depth: continue

      break

   if i == (max_attempt - 1):
      assert gMom.getHeight() <= max_depth
      return (gMom, gDad)
   else:
      nodeMom, nodeDad = momRandom, dadRandom

   nodeMom_parent = nodeMom.getParent()
   nodeDad_parent = nodeDad.getParent()

   # Sister
   if args["count"] >= 1:
      sister = gMom
      nodeDad.setParent(nodeMom_parent)

      if nodeMom_parent is None:
         sister.setRoot(nodeDad)
      else:
         nodeMom_parent.replaceChild(nodeMom, nodeDad)
      sister.processNodes()

      # self define, fixed ecoc tree
      _crossover_supplement_ecoc(sister, nodeDad)
      assert sister.getHeight() <= max_depth

   # Brother
   if args["count"] == 2:
      brother = gDad
      nodeMom.setParent(nodeDad_parent)

      if nodeDad_parent is None:
         brother.setRoot(nodeMom)
      else:
         nodeDad_parent.replaceChild(nodeDad, nodeMom)
      brother.processNodes()

      # self define, fixed ecoc tree
      _crossover_supplement_ecoc(brother, nodeMom)
      assert brother.getHeight() <= max_depth

   return (sister, brother)



