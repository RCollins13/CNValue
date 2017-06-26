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

#####Prepare annotation directory tree - anno sets
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

#####Prepare annotation directory tree - gene sets
mkdir ${WRKDIR}/analysis/geneSet_burden
while read pheno; do
  if [ -e ${WRKDIR}/analysis/geneSet_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/geneSet_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/geneSet_burden/${pheno}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}
    for filt in all; do
      if [ -e ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt} ]; then
        rm -rf ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}
      fi
      mkdir ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}
      for VF in E2 E3 E4 N1; do
        if [ -e ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF} ]; then
          rm -rf ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}
        fi
        mkdir ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}
        for context in exonic wholegene; do
          if [ -e ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/${context} ]; then
            rm -rf ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/${context}
          fi
          mkdir ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/${context}
        done
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Run _preliminary_ annotation set burden testing
while read pheno; do
  for CNV in CNV DEL DUP; do
    # for VF in E2 N1; do
    # for VF in E3 E4; do
    for VF in E2 E3 E4 N1; do
      # for filt in all noncoding; do
      for filt in haplosufficient; do
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

#####Run _secondary_ annotation set burden testing
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 N1; do
    # for VF in E3 E4; do
      # for filt in all noncoding; do
      for filt in haplosufficient; do
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

#####Collect annotation set burden test results
#Prepare directory
mkdir ${WRKDIR}/analysis/annoSet_burden/merged_results
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for filt in all haplosufficient noncoding; do
      for collection in effectSize pValue lowerCI upperCI zScore; do
        bsub -q short -sla miket_sc -u nobody -J ${CNV}_${VF}_${filt}_${collection} \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_annoSet_burdens.all_tests.sh \
        ${CNV} ${VF} ${filt} ${collection}"
      done
    done
  done
done
#Copy to plot data directory
if [ -e ${WRKDIR}/data/plot_data/annoSet_burden_results ]; then
  rm -rf ${WRKDIR}/data/plot_data/annoSet_burden_results
fi
cp -r ${WRKDIR}/analysis/annoSet_burden/merged_results \
${WRKDIR}/data/plot_data/annoSet_burden_results

#####Run gene set burden testing
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all; do
        #Exonic
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_geneSet_burdens_exonic \
        "${WRKDIR}/bin/rCNVmap/bin/geneSet_burdenTest_batch.sh -N 1000 \
          -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
          -p ${pheno}_${CNV}_${filt}_${VF}_exonic \
          -o ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/exonic/ \
          ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.list \
          ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf"
        #Whole-gene
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_geneSet_burdens_wholegene \
        "${WRKDIR}/bin/rCNVmap/bin/geneSet_burdenTest_batch.sh -N 1000 -W \
          -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
          -p ${pheno}_${CNV}_${filt}_${VF}_wholegene \
          -o ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/wholegene/ \
          ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.list \
          ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect gene set burden test results
#Prepare directory
mkdir ${WRKDIR}/analysis/geneSet_burden/merged_results
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for filt in all; do
      for context in exonic wholegene; do
        for collection in effectSize pValue lowerCI upperCI zScore; do
          bsub -q short -sla miket_sc -u nobody \
          -J ${CNV}_${VF}_${filt}_${context}_${collection} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneSet_burdens.all_tests.sh \
          ${CNV} ${VF} ${filt} ${context} ${collection}"
        done
      done
    done
  done
done
#Copy to plot data directory
if [ -e ${WRKDIR}/data/plot_data/geneSet_burden_results ]; then
  rm -rf ${WRKDIR}/data/plot_data/geneSet_burden_results
fi
cp -r ${WRKDIR}/analysis/geneSet_burden/merged_results \
${WRKDIR}/data/plot_data/geneSet_burden_results

#####Concatenate noncoding and coding burden test results into master matrixes
#Prepare directory
mkdir ${WRKDIR}/analysis/allSet_burden/merged_results
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for collection in effectSizes pvals lowerCI upperCI zScore; do
      for filt in all haplosufficient noncoding; do
        for dummy in 1; do
          #geneSet - filt=all for all noncoding filters
          #report wholegene for both deletions & duplications
          cat ${WRKDIR}/analysis/geneSet_burden/merged_results/${CNV}_${VF}_all_wholegene.${collection}.txt | \
          awk '{ print "GENESET."$0 }'
          #annoSet
          cat ${WRKDIR}/analysis/annoSet_burden/merged_results/${CNV}_${VF}_${filt}.${collection}.txt | \
          awk '{ print "ANNOSET."$0 }'
        done > ${WRKDIR}/analysis/allSet_burden/merged_results/${CNV}_${VF}_${filt}.${collection}.txt
      done
    done
  done
done
#Copy to plot data directory
if [ -e ${WRKDIR}/data/plot_data/allSet_burden_results ]; then
  rm -rf ${WRKDIR}/data/plot_data/allSet_burden_results
fi
cp -r ${WRKDIR}/analysis/allSet_burden/merged_results \
${WRKDIR}/data/plot_data/allSet_burden_results

#####Run exonic CNV burden testing for gene sets of constraint obs:exp deciles
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all; do
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_geneSet_burdens_exonic_constraint_deciles \
        "${WRKDIR}/bin/rCNVmap/bin/geneSet_burdenTest_batch.sh -f -N 1000 \
          -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
          -p ${pheno}_${CNV}_${filt}_${VF}_exonic \
          -o ${WRKDIR}/analysis/geneSet_burden/${pheno}/${CNV}/${filt}/${VF}/exonic/ \
          ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
          ${WRKDIR}/bin/rCNVmap/misc/LoF_constraint_decile_gene_sets.list \
          ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )







