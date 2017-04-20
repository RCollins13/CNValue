#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to parallelize annotation set CNV burden tests by phenotype-CNV combo

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
pheno=$1
CNV=$2

#####Run
for filt in all coding haplosufficient noncoding intergenic; do
  for VF in E2 E3 E4 N1; do
    echo -e "\n\n${pheno} ${CNV} ${filt} ${VF}\n\n"
    #Submit batch job
    ${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
      -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
      -p ${pheno}_${CNV}_${filt}_${VF} \
      -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list \
      ${WRKDIR}/data/misc/GRCh37_autosomes.genome
  done
done