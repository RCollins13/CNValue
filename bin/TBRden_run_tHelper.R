#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs binomial annotation helper R function

#Load library
library(TBRden)

#Read args
args <- commandArgs(trailingOnly=T)

#Run TBRden
tHelper(obs=args[1],perm=args[2],alt=args[3],OUTDIR=args[4],prefix=args[5])
