# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 19:51:46 2017

@author: Shone
"""

import Parsers, ParseColumn

from utils.stack import Stack


def _parseRunner(root, genids):
    
    ps = Stack()
    # ps.push(Parsers.ParserAHD())
    ps.push(ParseColumn.ParseColumn())
    # ps.push(Parsers.ParserFitness())
    ps.push(Parsers.ParserTestAcc())
    ps.push(Parsers.ParserTestFscore())
    ps.push(Parsers.ParserTrainAcc())
    ps.push(Parsers.ParserTrainFscore())
    
    while not ps.isEmpty():
        ps.pop().parse(root,genids)
