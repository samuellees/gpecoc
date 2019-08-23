#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2018/3/19 12:23
@author: Eric

To record good column of every individual in every generation.
It is used in local improvement.
# 
"""


class GeneBank():

    MAX_SPACE = 30
    MAX_ERROR = 0.5

    def __init__(self):
        self.genes = list()
        self.gene_num = 0

    def addgene(self, fcolumn, est_accuracy, class_accuracies, used_frequence):
        # duplicate check
        for gene in self.genes:
            if (fcolumn[1:] == gene[0][1:]).all() or (fcolumn[1:]*-1 == gene[0][1:]).all(): return
        # the estimator's max error rate should be lower than 0.5
        if 1 - est_accuracy < GeneBank.MAX_ERROR:
            self.genes.append((fcolumn, est_accuracy, class_accuracies, used_frequence))
            self.gene_num += 1
        # sort
        self.genes = sorted(self.genes, key=lambda gene: (gene[3], gene[1]), reverse=True)
        # Check if touch space limit
        if self.gene_num > GeneBank.MAX_SPACE:
            # delete estimator with a lowest accuracy
            self.genes = self.genes[:-1]
            self.gene_num -= 1

    def clear(self):
        self.genes = list()
        self.gene_num = 0

