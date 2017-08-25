#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects per-anno urCNV burden scores for a list of noncoding element BEDs

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reads arguments
pheno=$1
CNV=$2
VF=$3
filt=$4
list=$5

#####Iterates over elements list and collects data
while read class path; do
  if ! [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed.gz ] || \
     ! [ -s ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed.gz ] || \
       [ $( zcat ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed.gz | wc -l ) \
         -lt $( zcat ${path} | wc -l ) ]; then
    ${WRKDIR}/bin/rCNVmap/bin/gather_annoScore_data.sh -z \
    -p ${class} \
    -o ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
    ${path} \
    ${h37}
    if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed ]; then
      gzip -f ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScoreData.bed
    fi
  else
    echo "OUTPUT FOR ${pheno} ${CNV} ${VF} ${filt} vs ${path} ALREADY EXISTS. SKIPPING..."
  fi
done < ${list}