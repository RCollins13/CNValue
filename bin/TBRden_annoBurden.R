#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs functional annotation burden script

#Load library
library(TBRden)

#Read args
args <- commandArgs(trailingOnly=T)

#Run TBRden
annoBurden(anno=args[1],testBins=args[2],measure=args[3],
	OUTDIR=args[4],prefix=args[5],plot=T)
