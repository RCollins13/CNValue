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

#####Create directory if necessary
if ! [ -e ${WRKDIR}/data/plot_data/suppTables ]; then
  mkdir ${WRKDIR}/data/plot_data/suppTables
fi

#####Supp tables 1-2: significant large segments
VF=E2; filt=all
for CNV in DEL DUP; do
  #Print header
  for wrapper in 1; do
    echo -e "chr\tstart\tend\tCTRL_${CNV}_count"
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      echo -e "${pheno}_${CNV}_count\t${pheno}_${CNV}_OR\t${pheno}_${CNV}_p"
    done
  done | paste -s > ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt
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
  paste ${TMPDIR}/${CNV}_segments_stats.txt ${TMPDIR}/${CNV}_segments_CNVcounts.txt | \
  awk -v OFS="\t" '{ print $1, $2, $3, $14, $15, $4, $9, $16, $5, $10, $17, $6, $11, $18, $7, $12, $19, $8, $13 }' >> \
  ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt
done

