#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate organ-level master noncoding annotation tracks

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
tissue=$1

#####Run
while read anno; do
  echo ${anno}
  while read file; do
    cat ${file}
  done < <( awk -v tissue=${tissue} -v anno=${anno} \
  '{ if ($1==tissue && $2==anno) print $3 }' \
  ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list ) | \
  sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/${tissue}_MASTER.${anno}.elements.bed
done < <( awk -v tissue=${tissue} '{ if ($1==tissue) print $2 }' \
  ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
  sort | uniq )