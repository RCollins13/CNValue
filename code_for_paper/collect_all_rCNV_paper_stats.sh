#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master code log for all counts/stats in the paper

#Load parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Count total number of CNVs in final database
for filt in all noncoding; do
  echo ${filt}
  for group in GERM CNCR CTRL; do
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.${filt}.bed.gz | \
    fgrep -v "#" | wc -l
  done | awk '{ sum+=$1 }END{ print sum }'
done | paste - -
#Noncoding
for group in GERM CNCR CTRL; do
  zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.noncoding.bed.gz | \
  fgrep -v "#" | wc -l
done | awk '{ sum+=$1 }END{ print sum }'

#####Count number of coding & noncoding annotation classes
for class in gene_sets noncoding_annotations; do
  echo ${class}
  cat ${WRKDIR}/bin/rCNVmap/misc/master_${class}.list | wc -l
done | paste - -
#Per organ group
while read organ; do
  echo ${organ}
  for class in gene_sets noncoding_annotations; do
    grep -e "^${organ}" ${WRKDIR}/bin/rCNVmap/misc/master_${class}.list | wc -l
  done
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
  sort | uniq ) | paste - - - | awk -v OFS="\t" '{ print $1, $2, $3, $2+$3 }'

#####Get average number of CNVs per germline and cancer genome
#Germline
for group in CTRL GERM; do
  zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.all.bed.gz
done | fgrep -v "#" | wc -l | awk '{ print $1/102257 }'
#Cancer
zcat ${WRKDIR}/data/CNV/CNV_MASTER/CNCR/CNCR.CNV.E2.GRCh37.all.bed.gz | \
fgrep -v "#" | wc -l | awk '{ print $1/10844 }'

#####Get median CNV size by pheno
for group in CTRL GERM CNCR; do
  echo ${group}
  zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.all.bed.gz | \
  fgrep -v "#" | awk '{ print $3-$2 }' | sort -nk1,1 | \
  perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
done | paste - -
#All germline
for group in CTRL GERM; do
  zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.all.bed.gz | \
  fgrep -v "#" | awk '{ print $3-$2 }' 
done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
