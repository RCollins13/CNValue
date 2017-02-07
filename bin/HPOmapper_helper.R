#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Runs HPOmapper function via Rscript

#Load library
library(TBRden)

#Read args: 1) phenotypes, 2) HPO mappings, 3) output file
args <- commandArgs(trailingOnly=T)

#Run TBRden
HPOmapper(args[1],args[2],OUTFILE=args[3],gzip=T)