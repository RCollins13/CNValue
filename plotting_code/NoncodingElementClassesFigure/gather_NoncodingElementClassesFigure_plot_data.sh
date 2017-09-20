#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure for noncoding element class-based rCNV burden tests

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure ]; then
  rm -rf ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure
fi
mkdir ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure

#####Create matrix of fold-enrichments per annotation class for all phenotypes
VF=E2
filt=haplosufficient
if ! [ -e ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/burdenMatrices ]; then
  mkdir ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/burdenMatrices
fi
for CNV in DEL DUP; do
  cp ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.*.txt \
  ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/burdenMatrices/
done
