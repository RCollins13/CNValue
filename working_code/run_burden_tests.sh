#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master burden test code
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Prepare annotation directory tree
mkdir ${WRKDIR}/analysis/annoSet_burden
while read pheno; do
  if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}
    for filt in all coding haplosufficient noncoding intergenic; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt} ]; then
        rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}
      fi
      mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}
      for VF in E2 E3 E4 N1; do
        if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF} ]; then
          rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}
        fi
        mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Run _preliminary_ annotation set burden testing
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 N1; do
      for filt in all noncoding; do
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_annoSet_burdens \
        "${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
          -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
          -p ${pheno}_${CNV}_${filt}_${VF} \
          -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
          ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prelim_subset.list \
          ${WRKDIR}/data/misc/GRCh37_autosomes.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect _preliminary annotation set burden test results
#Prepare directory
mkdir ${WRKDIR}/analysis/annoSet_burden/merged_results
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 N1; do
    for filt in all noncoding; do
      for collection in effectSize pValue lowerCI upperCI zScore; do
        bsub -q short -sla miket_sc -u nobody -J ${CNV}_${VF}_${filt}_${collection} \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_annoSet_burdens.preliminary_tests.sh \
        ${CNV} ${VF} ${filt} ${collection}"
      done
    done
  done
done
#Copy to plot data directory
cp -r ${WRKDIR}/analysis/annoSet_burden/merged_results \
${WRKDIR}/data/plot_data/annoSet_burden_results

#####Run _secondary_ annotation set burden testing
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 N1; do
      for filt in all noncoding; do
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_annoSet_burdens_secondary \
        "${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
          -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
          -p ${pheno}_${CNV}_${filt}_${VF} \
          -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
          ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.secondary_subset.list \
          ${WRKDIR}/data/misc/GRCh37_autosomes.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

# #####Run _full_ annotation set burden testing
# #iterate & prep directories
# while read pheno; do
#   for CNV in CNV DEL DUP; do
#     #CODE:
#     # for filt in all coding haplosufficient noncoding intergenic; do
#     #   for VF in E2 E3 E4 N1; do
#     #     #Submit batch job
#     #     ${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
#     #       -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
#     #       -p ${pheno}_${CNV}_${filt}_${VF} \
#     #       -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
#     #       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
#     #       ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
#     #       ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list \
#     #       ${WRKDIR}/data/misc/GRCh37_autosomes.genome
#     #   done
#     # done
#     #PARALLELIZE:
#     bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_annoSet_burdens \
#     "${WRKDIR}/bin/rCNVmap/analysis_scripts/parallelize_batched_annoSet_burdens.sh \
#     ${pheno} ${CNV}"
#   done
# done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
#           cut -f1 | fgrep -v CTRL )




