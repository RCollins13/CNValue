#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs core TBRden function

#Load library
library(TBRden)

#Read args
args <- commandArgs(trailingOnly=T)

#Run TBRden
TBRden(controls=as.character(args[1]),
	   cases=as.character(args[2]),
	   OUTDIR=as.character(args[3]),
	   prefix=as.character(args[4]),
	   adjusted=as.numeric(args[5]),
	   manColor=as.character(args[6]))
