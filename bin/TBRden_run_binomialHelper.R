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
binomialHelper(obs=args[1],n=args[2],perm=args[3],
	alt=args[4],OUTDIR=args[5],prefix=args[6])
