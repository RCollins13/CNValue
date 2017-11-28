#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for figure estimating risk across phenotypes & classes

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/riskEstimatesFigure ]; then
  rm -rf ${WRKDIR}/data/plot_data/riskEstimatesFigure
fi
mkdir ${WRKDIR}/data/plot_data/riskEstimatesFigure

#####Create temporary file with estimates of regulatory elements to consider
for CNV in DEL DUP; do
  cut -f5 ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
  sed -e 's/\;/\n/g' -e 's/\_/\t/g' | sort -Vk1,1 -k2,2n -k3,3n | \
  grep -e '^[0-9]' | bedtools merge -i - > \
  ${TMPDIR}/noncoding_${CNV}.sites.bed
done

#####Risk estimates per phenotype per class
#Large loci
while read pheno ncase; do
  for wrapper in 1; do
    echo "${pheno}"
    for CNV in DEL DUP; do
      #CNV incidence
      sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | \
      cut -f1-3 | bedtools coverage -a - \
      -b <( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E2.GRCh37.all.bed.gz | \
            awk -v OFS="\t" '{ if ($3-$2>=200000) print $1, $2, $3 }' ) | \
      awk '{ if ($(NF-2)>=200000) print $4 }' | wc -l | \
      awk -v ncase=${ncase} '{ print $1/ncase }'
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1,8 | grep -e 'CTRL\|GERM\|NEURO\|NDD\|PSYCH\|SOMA' ) > \
${WRKDIR}/data/plot_data/riskEstimatesFigure/largeSegment.incidence.txt
#Genes
while read pheno ncase; do
  for wrapper in 1; do
    echo "${pheno}"
    for CNV in DEL DUP; do
      #CNV incidence
      sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | \
      cut -f1 | sed 's/\-/_/g' | fgrep -wf - \
      <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed ) | \
      bedtools intersect -u -b - \
      -a ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.all_largeSegmentExcluded.bed.gz | \
      wc -l | awk -v ncase=${ncase} '{ print $1/ncase }'
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1,8 | grep -e 'CTRL\|GERM\|NEURO\|NDD\|PSYCH\|SOMA' ) > \
${WRKDIR}/data/plot_data/riskEstimatesFigure/gene.incidence.txt
#Regulatory blocks
while read pheno ncase; do
  for wrapper in 1; do
    echo "${pheno}"
    for CNV in DEL DUP; do
      #Set appropriate CNV file
      if [ ${CNV} == "DEL" ]; then
        CNVfile=${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient_largeSegmentExcluded_sigGeneExcluded.bed.gz
      elif [ ${CNV} == "DUP" ]; then
        CNVfile=${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.noncoding_largeSegmentExcluded.bed.gz
      fi
      #CNV incidence
      bedtools intersect -u -a ${CNVfile} -b ${TMPDIR}/noncoding_${CNV}.sites.bed | \
      wc -l | awk -v ncase=${ncase} '{ print $1/ncase }'
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1,8 | grep -e 'CTRL\|GERM\|NEURO\|NDD\|PSYCH\|SOMA' ) > \
${WRKDIR}/data/plot_data/riskEstimatesFigure/regBlock.incidence.txt


