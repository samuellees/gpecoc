import random

import Util
from GenomeBase import GenomeBase, GTreeBase, GTreeNodeBase


nodeType = {"TERMINAL" : 0, "NONTERMINAL": 1}

class GTreeNodeGP(GTreeNodeBase):
   """ The GTreeNodeGP Class - The Genetic Programming Node representation
   
   Inheritance diagram for :class:`GTree.GTreeNodeGP`:

   .. inheritance-diagram:: GTree.GTreeNodeGP

   :param data: the node data
   :param type: the node type
   :param parent: the node parent
   
   """
   def __init__(self, data, node_type=0, parent=None):
      GTreeNodeBase.__init__(self, parent)
      self.node_type = node_type
      self.node_data = data

   def __repr__(self):
      str_repr  = GTreeNodeBase.__repr__(self)
      str_repr += " - [%s]" % self.node_data
      return str_repr     

   def compare(self, other):
      """ Compare this node with other 
      
      :param other: the other GTreeNodeGP
      """
      if not isinstance(other, GTreeNodeGP):
         Util.raiseException("The other node used to compare is not a GTreeNodeGP class", TypeError)

      if other.node_type == self.node_type:
         if other.node_data == self.node_data:
            return 0
      return -1

   def setData(self, data):
      """Sets the node internal data
      
      :param data: the internal data
      """
      self.node_data = data

   def getData(self):
      """Gets the node internal data
      
      :rtype: the internal data
      """
      return self.node_data

   def setType(self, node_type):
      """Sets the node type 
      
      :param node_type: the node type is type of Consts.nodeType
      """
      self.node_type = node_type

   def getType(self):
      """Get the node type 
      
      :rtype: the node type is type of Consts.nodeType
      """
      return self.node_type

   def newNode(self, data):
      """Creates a new node and adds this
      node as children of current node

      :param data: the internal node data
      """
      node = GTreeNodeGP(data, self)
      self.addChild(node)
      return node

   def swapNodeData(self, node):
      """Swaps the node data and type with another node

      :param node: the node
      """
      tmp_data = self.node_data
      tmp_type = self.node_type
      self.setData(node.getData())
      self.setType(node.getType())
      node.setData(tmp_data)
      node.setType(tmp_type)

   def copy(self, g):
      """ Copy the contents to the destination g
      
      :param g: the GTreeNodeGP genome destination
      """
      GTreeNodeBase.copy(self, g)
      g.node_data = self.node_data
      g.node_type = self.node_type

   def clone(self):
      """ Return a new copy of the node

      :rtype: the new GTreeNodeGP instance
      """
      newcopy = GTreeNodeGP(None)
      self.copy(newcopy)
      return newcopy



def checkTerminal(terminal):
   """ Do some check on the terminal, to evaluate ephemeral constants

   :param terminal: the terminal string
   """
   if terminal.startswith("ephemeral:"):
      splited = terminal.split(":")
      ephemeral_constant = eval(splited[1])
      return str(ephemeral_constant)
   else:
      return terminal

def buildGTreeGPGrow(ga_engine, depth, max_depth):
   """ Creates a new random GTreeGP root node with subtrees using
   the "Grow" method.
   
   :param ga_engine: the GA Core
   :param depth: the initial depth
   :max_depth: the maximum depth of the tree
   :rtype: the root node
   """

   gp_terminals = ga_engine.getParam("gp_terminals")
   assert gp_terminals is not None

   gp_function_set = ga_engine.getParam("gp_function_set")
   assert gp_function_set is not None

   if depth == max_depth:
      random_terminal = checkTerminal(random.choice(gp_terminals))
      n = GTreeNodeGP(random_terminal, nodeType["TERMINAL"])
      return n
   else:
      # Do not generate degenerative trees 
      if depth == 0:
         random_node = random.choice(gp_function_set.keys())
      else:
         fchoice = random.choice([gp_function_set.keys(), gp_terminals])
         random_node = random.choice(fchoice)

      if random_node in gp_terminals:
         n = GTreeNodeGP(checkTerminal(random_node), nodeType["TERMINAL"])
      else:
         n = GTreeNodeGP(random_node, nodeType["NONTERMINAL"])

   if n.getType() == nodeType["NONTERMINAL"]:
      for i in xrange(gp_function_set[n.getData()]):
         child = buildGTreeGPGrow(ga_engine, depth+1, max_depth)
         child.setParent(n)
         n.addChild(child)

   return n

def buildGTreeGPFull(ga_engine, depth, max_depth):
   """ Creates a new random GTreeGP root node with subtrees using
   the "Full" method.
   
   :param ga_engine: the GA Core
   :param depth: the initial depth
   :max_depth: the maximum depth of the tree
   :rtype: the root node
   """
   gp_terminals = ga_engine.getParam("gp_terminals")
   assert gp_terminals is not None

   gp_function_set = ga_engine.getParam("gp_function_set")
   assert gp_function_set is not None

   if depth == max_depth:
      random_terminal = checkTerminal(random.choice(gp_terminals))
      n = GTreeNodeGP(random_terminal, nodeType["TERMINAL"])
      return n
   else:
      random_oper = random.choice(gp_function_set.keys())
      n = GTreeNodeGP(random_oper, nodeType["NONTERMINAL"])

   if n.getType() == nodeType["NONTERMINAL"]:
      for i in xrange(gp_function_set[n.getData()]):
         child = buildGTreeGPFull(ga_engine, depth+1, max_depth)
         child.setParent(n)
         n.addChild(child)

   return n

