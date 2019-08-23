"""

:mod:`Consts` -- constants module
============================================================================

Pyevolve have defaults in all genetic operators, settings and etc, this is an issue to helps the user in the API use and minimize the source code needed to make simple things. In the module :mod:`Consts`, you will find those defaults settings. You are encouraged to see the constants, but not to change directly on the module, there are methods for this.

General constants
----------------------------------------------------------------------------

.. attribute:: CDefPythonRequire
  
   The mininum version required to run Pyevolve.

.. attribute:: CDefLogFile
   
   The default log filename.

.. attribute:: CDefLogLevel

   Default log level.

.. attribute:: sortType
   
   Sort type, raw or scaled.

   Example:
      >>> sort_type = Consts.sortType["raw"]
      >>> sort_type = Consts.sortType["scaled"]

.. attribute:: minimaxType

   The Min/Max type, maximize or minimize the evaluation function.

   Example:
      >>> minmax = Consts.minimaxType["minimize"]
      >>> minmax = Consts.minimaxType["maximize"]
  
.. attribute:: CDefESCKey

   The ESC key ASCII code. Used to start Interactive Mode.

.. attribute:: CDefRangeMin

   Minimum range. This constant is used as integer and real max/min.

.. attribute:: CDefRangeMax

   Maximum range. This constant is used as integer and real max/min.

.. attribute:: CDefBroadcastAddress
   
   The broadcast address for UDP, 255.255.255.255

.. attribute:: CDefImportList
   
   The import list and messages

.. attribute:: nodeType

   The genetic programming node types, can be "TERMINAL":0 or "NONTERMINAL":1

.. attribute:: CDefGPGenomes

   The classes which are used in Genetic Programming, used to detected the
   correct mode when starting the evolution

Selection methods constants (:mod:`Selectors`)
----------------------------------------------------------------------------

.. attribute:: CDefTournamentPoolSize

   The default pool size for the Tournament Selector (:func:`Selectors.GTournamentSelector`).

Scaling scheme constants (:mod:`Scaling`)
----------------------------------------------------------------------------

.. attribute:: CDefScaleLinearMultiplier

   The multiplier of the Linear (:func:`Scaling.LinearScaling`) scaling scheme.

.. attribute:: CDefScaleSigmaTruncMultiplier

   The default Sigma Truncation (:func:`Scaling.SigmaTruncScaling`) scaling scheme.

.. attribute:: CDefScalePowerLawFactor

   The default Power Law (:func:`Scaling.PowerLawScaling`) scaling scheme factor.

.. attribute:: CDefScaleBoltzMinTemp

   The default mininum temperature of the (:func:`Scaling.BoltzmannScaling`) scaling scheme factor.

.. attribute:: CDefScaleBoltzFactor

   The default Boltzmann Factor of (:func:`Scaling.BoltzmannScaling`) scaling scheme factor.
   This is the factor that the temperature will be subtracted.

.. attribute:: CDefScaleBoltzStart

   The default Boltzmann start temperature (:func:`Scaling.BoltzmannScaling`).
   If you don't set the start temperature parameter, this will be the default initial
   temperature for the Boltzmann scaling scheme.

Population constants (:class:`GPopulation.GPopulation`)
----------------------------------------------------------------------------
   
.. attribute:: CDefPopSortType
   
   Default sort type parameter.

.. attribute:: CDefPopMinimax

   Default min/max parameter.

.. attribute:: CDefPopScale

   Default scaling scheme.

Tree chromosome constants (:class:`GTree.GTree`)
----------------------------------------------------------------------------

.. attribute:: CDefGTreeInit

   Default initializator of the tree chromosome.
   
.. attribute:: CDefGGTreeMutator

   Default mutator of the tree chromosome.
   
.. attribute:: CDefGTreeCrossover

   Default crossover of the tree chromosome.

GA Engine constants (:class:`GSimpleGA.GSimpleGA`)
----------------------------------------------------------------------------

.. attribute:: CDefGAGenerations

   Default number of generations.

.. attribute:: CDefGAMutationRate

   Default mutation rate.

.. attribute:: CDefGACrossoverRate

   Default crossover rate.

.. attribute:: CDefGAPopulationSize

   Default population size.

.. attribute:: CDefGASelector

   Default selector method.


"""
import logging

import Crossovers
import Mutators
import Scaling
import Selectors
import Initializators

# module Consts
CDefLogFile = "pyevolve.log"
CDefLevel = logging.DEBUG

# Types of sort
# - raw: uses the "score" attribute
# - scaled: uses the "fitness" attribute
sortType = { 
   "raw"    : 0,
   "scaled" : 1
}

# Optimization type
# - Minimize or Maximize the Evaluator Function
minimaxType = {
    "minimize": 0,
    "maximize": 1
}

CDefESCKey = 27




####################
# Defaults section #
####################

# - Tournament selector
CDefTournamentPoolSize = 2

# - Scale methods defaults
CDefScaleLinearMultiplier     = 1.2
CDefScaleSigmaTruncMultiplier = 2.0
CDefScalePowerLawFactor       = 1.0005
CDefScaleBoltzMinTemp         = 1.0
CDefScaleBoltzFactor          = 0.05
# 40 temp. = 500 generations
CDefScaleBoltzStart           = 40.0

# - Population Defaults
CDefPopSortType               = sortType["raw"]
CDefPopMinimax                = minimaxType["maximize"]
CDefPopScale                  = Scaling.LinearScaling

# - GA Engine defaults
CDefGAGenerations    = 100
CDefGAMutationRate   = 0.02
CDefGACrossoverRate  = 0.9
CDefGAPopulationSize = 80




#CDefGASelector       = Selectors.GRankSelector
#CDefGASelector       = Selectors.GRouletteWheel
CDefGASelector       = Selectors.GTournamentSelectorAlternative
#CDefGASelector       = Selectors.GTournamentSelector




CDefGAElitismReplacement = 1

# - This is general used by integer/real ranges defaults
CDefRangeMin = 0
CDefRangeMax = 100

# - GTreeGP defaults
CDefGTreeGPInit      = Initializators.GTreeGPInitializator


#CDefGGTreeGPMutator  = Mutators.GTreeGPMutatorSubtree
# CDefGGTreeGPMutator  = Mutators.DIYGTreeGPMutatorSubtree
CDefGGTreeGPMutator  = Mutators.selfDefined_GTreeGPMutatorSubtree


# CDefGTreeGPCrossover = Crossovers.GTreeGPCrossoverSinglePoint
CDefGTreeGPCrossover = Crossovers.NewGTreeGPCrossover_ECOC

# Util Consts
CDefBroadcastAddress = "255.255.255.255"
nodeType = {"TERMINAL" : 0, "NONTERMINAL": 1}



