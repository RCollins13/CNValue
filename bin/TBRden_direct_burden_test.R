#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs TBRden direct test function

#Load library
library(TBRden)

#Read args
args <- commandArgs(trailingOnly=T)

#Order:
#1: observed
#2: permuted
#3: tail
#4: manColor
#5: OUTDIR
#6: prefix

#Run TBRden
directTest(observed=args[1],
           permuted=args[2],
           tail=args[3],
           Bonferroni=T,
           QQ=T,
           manhattan=T,
           manMaxY=50,
           manColor=args[4],
           return=F,
           OUTDIR=args[5],
           prefix=args[6],
           gzip=T)