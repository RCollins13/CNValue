#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to compile all supplementary tables for rCNV paper

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Supp tables 1-2: significant large segments
VF=E2; filt=all
for CNV in DEL DUP; do
  bedtools intersect -wa -wb \
  -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR.bed | \
        cut -f4-6,10,16 ) | \
  bedtools intersect -wa -wb -a - \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal.bed | \
        cut -f4-6,10,16 )
done

