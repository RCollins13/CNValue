#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Gets odds ratios for a list of regulatory blocks

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
VF=$1
filt=$2
CNV=$3

#####Run
while read chr start end blockID elements eIDs; do
  for dummy in 1; do
    echo -e "chr${chr}\n${start}\n${end}\n${blockID}"
    #Iterate over phenos
    while read pheno nCASE; do
      #Determine CNV files
      if [ ${CNV} == "DEL" ]; then
        CTRL_CNVs=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.haplosufficient_largeSegmentExcluded_sigGeneExcluded.bed.gz
        CASE_CNVs=${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.haplosufficient_largeSegmentExcluded_sigGeneExcluded.bed.gz
      elif [ ${CNV} == "DUP" ]; then
        CTRL_CNVs=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.noncoding_largeSegmentExcluded.bed.gz
        CASE_CNVs=${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.noncoding_largeSegmentExcluded.bed.gz
      fi
      #Get counts of case/control CNVs
      caseCNV=$( echo ${elements} | sed -e 's/\;/\n/g' -e 's/_/\t/g' | \
      bedtools intersect -u -b - -a ${CASE_CNVs} | cut -f4 | sort | uniq | wc -l )
      controlCNV=$( echo ${elements} | sed -e 's/\;/\n/g' -e 's/_/\t/g' | \
      bedtools intersect -u -b - -a ${CTRL_CNVs} | cut -f4 | sort | uniq | wc -l )
      caseNoCNV=$( echo "${nCASE}-${caseCNV}" | bc )
      controlNoCNV=$( echo "38628-${controlCNV}" | bc )
      if [ ${pheno} == "GERM" ]; then
        echo "${controlCNV}"
      fi
      echo -e "${caseCNV}"
      #Calcluate odds ratio
      unset R_HOME
      Rscript -e "cat(paste(round(((${caseCNV}+(${nCASE}/38628))/(${controlCNV}+1))/((${caseNoCNV}+(${nCASE}/38628))/(${controlNoCNV}+1)),4)),\"\n\",sep=\"\")"
      #Calcluate p-value
      Rscript -e "cat(paste(format(fisher.test(matrix(c(${controlNoCNV},${caseNoCNV},${controlCNV},${caseCNV}),byrow=T,nrow=2),alternative=\"greater\")\$p.value,scientific=T)),\"\n\",sep=\"\")"
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      cut -f1,8 | grep -e 'GERM\|NEURO\|NDD\|PSYCH\|SOMA' )
  done | paste -s
done < ${WRKDIR}/analysis/perAnno_burden/signifRegulatoryBlocks.final.bed > \
${WRKDIR}/analysis/perAnno_burden/signifRegulatoryBlocks.${CNV}_ORs.bed
