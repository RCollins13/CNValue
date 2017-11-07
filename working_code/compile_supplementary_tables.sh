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
  #Get ORs and p-values
  bedtools intersect -wa -wb \
  -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR.bed ) | \
  bedtools intersect -wa -wb -a - \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal.bed ) | \
  cut --complement -f4-6,12-14 > \
  ${TMPDIR}/${CNV}_segments_stats.txt
  #Get count of CNVs per pheno
  while read chr start end; do
    for pheno in CTRL GERM NEURO NDD PSYCH SOMA; do
      bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
      cut -f4- | bedtools coverage -b - -a <( echo -e "${chr}\t${start}\t${end}" ) | \
      awk -v overlap=${overlap} '{ if ($(NF-2)>=overlap) print $4 }' | sort | uniq | wc -l
    done | paste -s
  done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed > \
  ${TMPDIR}/${CNV}_segments_CNVcounts.txt
  #Paste the two tables together
  paste ${TMPDIR}/${CNV}_segments_stats.txt ${TMPDIR}/${CNV}_segments_CNVcounts.txt
done

