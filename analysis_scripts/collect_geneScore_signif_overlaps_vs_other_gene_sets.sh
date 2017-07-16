#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects gene count overlap of autosomal protein-coding genes from unique signif
# per-gene tests (geneScore) vs a supplied set of comparison gene lists

# Automatically runs for CTRL, GERM, NEURO, NDD, PSYCH, SOMA, CNCR (in order)

#Read arguments
CNV=$1
VF=$2
context=$3
signif=$4
complist=$5 #two columns: gene list name & full path. Tab-delimmed

