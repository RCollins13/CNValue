#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run all rCNV burden scoring per annotation

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directories if exists
if [ -e ${WRKDIR}/data/perAnno_burden ]; then
  rm -rf ${WRKDIR}/data/perAnno_burden
fi
mkdir ${WRKDIR}/data/perAnno_burden
if [ -e ${WRKDIR}/analysis/perAnno_burden ]; then
  rm -rf ${WRKDIR}/analysis/perAnno_burden
fi
mkdir ${WRKDIR}/analysis/perAnno_burden

#####Create subdirectories
while read pheno; do
  if [ -e ${WRKDIR}/data/perAnno_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perAnno_burden/${pheno}
  if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Reorder noncoding annotations based on number of elements in group
while read class path; do
  nElements=$( cat ${path} | wc -l )
  echo -e "${class}\t${path}\t${nElements}"
done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list | \
sort -nk3,3 | cut -f1-2 > \
${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list

#####Submit burden data collection for all phenotypes
while read pheno; do
  for CNV in CNV DEL DUP; do
    # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV} ]; then
    #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}
    # fi
    # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}
    for VF in E4; do
      # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF} ]; then
      #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}
      # fi
      # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}
      for filt in haplosufficient noncoding; do
        # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt} ]; then
        #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        # fi
        # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perAnno_burden_dataCollection_${filt} \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/gather_annoScore_data_batchMode.sh -F \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Submit annoScore model for all phenotypes
while read pheno; do
  #Get number of subjects in group
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
           ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}
    for VF in E4; do
      if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF} ]; then
        rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}
      fi
      mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}
      for filt in haplosufficient noncoding; do
        if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt} ]; then
          rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        fi
        mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        #Launch model in batch mode across all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_annoScoreModel \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/run_annoScore_model_batchMode.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect significant elements per track per phenotype
while read pheno; do
  #Make output directories (if necessary)
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}
  fi
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL
  fi
  for CNV in CNV DEL DUP; do
    for VF in E4; do
      for filt in haplosufficient noncoding; do
        #Launch collection script for all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_annoScoreModel \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_perClass_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect significant elements across all tracks per phenotype
#Make lists of tissue-defined elements
while read tissue; do
  fgrep ${tissue} ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list > \
  ${WRKDIR}/lists/all_${tissue}_genome_annotations.list
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
          sort | uniq )
#Make list of tissue-agnostic annotations
cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
sort | uniq | fgrep -vf - \
${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list > \
${WRKDIR}/lists/tissue_agnostic_genome_annotations.list
#Iterate over phenotypes
while read pheno; do
  #Make output directories (if necessary)
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged
  fi
  for CNV in CNV DEL DUP; do
    for VF in E4; do
      for filt in haplosufficient noncoding; do
        #Launch merge & filtering script across all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_all \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list \
        bonf \
        ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.all_classes.bed"
        #Launch merge & filtering script for tissue-agnostic annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_tissueAgnostic \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/lists/tissue_agnostic_genome_annotations.list \
        bonf \
        ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.tissue_agnostic.bed"
        #Iterate over tissues and launch merge & filtering script for tissue-dependent annotations
        while read tissue; do
          bsub -q short -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_${tissue} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
          ${pheno} ${CNV} ${VF} ${filt} \
          ${WRKDIR}/lists/all_${tissue}_genome_annotations.list \
          bonf \
          ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.${tissue}.bed"
        done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq )
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )












