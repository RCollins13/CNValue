#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Runs per-anno urCNV burden model for a list of noncoding element BEDs

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reads arguments
pheno=$1
CNV=$2
VF=$3
filt=$4
list=$5

#####Set number of subjects per phenotype group
nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
         ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
nCTRL=38628

#####Iterates over elements list and runs model
while read class path; do
  annoScoreData=${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed.gz
  if [ -e ${annoScoreData} ]; then
    OUTFILE=${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScore_stats.bed
    #Runs model from gzipped file
    ${WRKDIR}/bin/rCNVmap/bin/run_annoScore_model.R \
    -o ${OUTFILE} \
    ${annoScoreData} \
    ${nCTRL} \
    ${nCASE}
    #Compress output
    if [ -e ${OUTFILE} ]; then
      gzip -f ${OUTFILE}
    fi
  else
    echo "OUTPUT FOR ${pheno} ${CNV} ${VF} ${filt} vs ${path} IS MISSING. SKIPPING..."
  fi
done < ${list}