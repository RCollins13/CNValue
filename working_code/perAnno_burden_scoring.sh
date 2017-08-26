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
        bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perAnno_burden_dataCollection_${filt} \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/gather_annoScore_data_batchMode.sh \
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
        #ADD COMMAND HERE
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

