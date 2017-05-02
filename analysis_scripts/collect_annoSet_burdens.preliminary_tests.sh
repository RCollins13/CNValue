#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to collect annotation set CNV burden test results

#Master burden test code
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
CNV=$1
VF=$2
filt=$3
collection=$4

#####Run
if [ ${collection} == "effectSize" ]; then
  while read anno annopath; do
    echo "${anno}"
    for pheno in GERM UNK NEURO NDD DD PSYCH SCZ ASD SEIZ HYPO BEHAV ID \
    SOMA HEAD GRO CARD SKEL DRU MSC EE INT EMI CNCR CGEN CSKN CGST \
    CRNL CBRN CLNG CBST CEND CHNK CLIV CMSK CBLD; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt ]; then
        nfold=$( fgrep -v "#" \
          ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt | \
          awk -v OFS="\t" '{ print $11 }' )
      else
        nfold=NA
      fi
      if [ -z ${nfold} ]; then
        nfold=NA
      fi
      echo ${nfold}
    done | paste -s
  done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prelim_subset.sorted.list | \
  paste - - > ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.effectSizes.txt
elif [ ${collection} == "pValue" ]; then
  while read anno annopath; do
    echo "${anno}"
    for pheno in GERM UNK NEURO NDD DD PSYCH SCZ ASD SEIZ HYPO BEHAV ID \
    SOMA HEAD GRO CARD SKEL DRU MSC EE INT EMI CNCR CGEN CSKN CGST \
    CRNL CBRN CLNG CBST CEND CHNK CLIV CMSK CBLD; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt ]; then
        p=$( fgrep -v "#" \
          ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt | \
          awk -v OFS="\t" '{ print $NF }' )
      else
        p=NA
      fi
      if [ -z ${p} ]; then
        p=NA
      fi
      echo ${p}
    done | paste -s
  done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prelim_subset.sorted.list | \
  paste - - > ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.pvals.txt
elif [ ${collection} == "lowerCI" ]; then
  while read anno annopath; do
    echo "${anno}"
    for pheno in GERM UNK NEURO NDD DD PSYCH SCZ ASD SEIZ HYPO BEHAV ID \
    SOMA HEAD GRO CARD SKEL DRU MSC EE INT EMI CNCR CGEN CSKN CGST \
    CRNL CBRN CLNG CBST CEND CHNK CLIV CMSK CBLD; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt ]; then
        CI=$( fgrep -v "#" \
          ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt | \
          awk -v OFS="\t" '{ print $12 }' )
      else
        CI=NA
      fi
      if [ -z ${CI} ]; then
        CI=NA
      fi
      echo ${CI}
    done | paste -s
  done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prelim_subset.sorted.list | \
  paste - - > ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.lowerCI.txt
elif [ ${collection} == "upperCI" ]; then
  while read anno annopath; do
    echo "${anno}"
    for pheno in GERM UNK NEURO NDD DD PSYCH SCZ ASD SEIZ HYPO BEHAV ID \
    SOMA HEAD GRO CARD SKEL DRU MSC EE INT EMI CNCR CGEN CSKN CGST \
    CRNL CBRN CLNG CBST CEND CHNK CLIV CMSK CBLD; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt ]; then
        CI=$( fgrep -v "#" \
          ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/${pheno}_${CNV}_${filt}_${VF}.${anno}.CNV_burden_results.txt | \
          awk -v OFS="\t" '{ print $13 }' )
      else
        CI=NA
      fi
      if [ -z ${CI} ]; then
        CI=NA
      fi
      echo ${CI}
    done | paste -s
  done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prelim_subset.sorted.list | \
  paste - - > ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.upperCI.txt
fi
