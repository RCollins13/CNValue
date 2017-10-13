#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in SMARCA2 figure

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/SMARCA2 ]; then
  rm -rf ${WRKDIR}/data/plot_data/SMARCA2
fi
mkdir ${WRKDIR}/data/plot_data/SMARCA2

#####Gather coordinates of NDD and CTRL CNVs in plotting window 
for pheno in NDD CTRL; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E4; do
      for filt in all coding haplosufficient; do
        zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
        fgrep -v "#"  | bedtools intersect -wa -a - -b <( echo -e "9\t1250000\t3250000" ) > \
        ${WRKDIR}/data/plot_data/SMARCA2/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed
        gzip -f ${WRKDIR}/data/plot_data/SMARCA2/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed
      done
    done
  done
done

#####Run CNV pileups in 100bp windows across entire plotting window
step=100
for pheno in NDD CTRL; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E4; do
      for filt in all coding haplosufficient; do
        paste <( seq 1250000 ${step} 3249900 ) <( seq 1250100 ${step} 3250000 ) | \
        awk -v OFS="\t" '{ print "9", $1, $2 }' | bedtools intersect -c -wa -a - \
        -b  ${WRKDIR}/data/plot_data/SMARCA2/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz > \
        ${WRKDIR}/data/plot_data/SMARCA2/${pheno}.${CNV}.${VF}.GRCh37.${filt}.100bp_pileups.bed
        gzip -f ${WRKDIR}/data/plot_data/SMARCA2/${pheno}.${CNV}.${VF}.GRCh37.${filt}.100bp_pileups.bed
      done
    done
  done
done
