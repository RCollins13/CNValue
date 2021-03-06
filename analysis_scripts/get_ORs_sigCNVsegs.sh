#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Gets odds ratios for a list of loci from the CNV hotspot analysis

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
VF=$1
filt=$2
CNV=$3
overlap=200000

#####Run
while read chr start end; do
  for dummy in 1; do
    echo -e "chr${chr}\n${start}\n${end}"
    #iterate over phenos
    while read pheno nCASE; do
      #Get counts of case/control CNVs
      caseCNV=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
      cut -f4- | bedtools coverage -b - -a <( echo -e "${chr}\t${start}\t${end}" ) | \
      awk -v overlap=${overlap} '{ if ($(NF-2)>=overlap) print $4 }' | sort | uniq | wc -l )
      controlCNV=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
      cut -f4- | bedtools coverage -b - -a <( echo -e "${chr}\t${start}\t${end}" ) | \
      awk -v overlap=${overlap} '{ if ($(NF-2)>=overlap) print $4 }' | sort | uniq | wc -l )
      caseNoCNV=$( echo "${nCASE}-${caseCNV}" | bc )
      controlNoCNV=$( echo "38628-${controlCNV}" | bc )
      #Calcluate odds ratio
      unset R_HOME
      Rscript -e "cat(paste(round(((${caseCNV}+(${nCASE}/38628))/(${controlCNV}+1))/((${caseNoCNV}+(${nCASE}/38628))/(${controlNoCNV}+1)),4)),\"\n\",sep=\"\")"
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      cut -f1,8 | grep -e 'GERM\|NEURO\|NDD\|PSYCH\|SOMA' )
  done | paste -s
done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed > \
${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR.bed
