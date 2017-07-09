#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run all rCNV burden scoring per gene

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/perGene_burden
fi
mkdir ${WRKDIR}/data/perGene_burden

#####Create subdirectories
while read pheno; do
  if [ -e ${WRKDIR}/data/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perGene_burden/${pheno}
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Submit burden data collection for all phenotypes
#Germline
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      #Exonic
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_exonic \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_exonic.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
      #Wholegene
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_wholegene \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh -W \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_wholegene.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Copy test datasets to plotting data directory
if [ -e ${WRKDIR}/data/plot_data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/plot_data/perGene_burden
fi
mkdir ${WRKDIR}/data/plot_data/perGene_burden
for pheno in GERM CNCR; do
  for CNV in CNV DEL DUP; do
    for VF in E2 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          cp ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          ${WRKDIR}/data/plot_data/perGene_burden/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt
        fi
      done
    done
  done
done

#####Copy all datasets to plotting data directory
if [ -e ${WRKDIR}/data/plot_data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/plot_data/perGene_burden
fi
mkdir ${WRKDIR}/data/plot_data/perGene_burden
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          cp ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          ${WRKDIR}/data/plot_data/perGene_burden/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Run geneScore model for all germline phenotype groups
#Create master output directory
if [ -e ${WRKDIR}/analysis/perGene_burden ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden
fi
mkdir ${WRKDIR}/analysis/perGene_burden
#Run model
while read pheno; do
  #Create output directory for results
  if [ -e ${WRKDIR}/data/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perGene_burden/${pheno}
  #Get number of subjects in group
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
           ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${context}_geneScoreModel \
          "${WRKDIR}/bin/rCNVmap/bin/run_geneScore_model.R \
          -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
          ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          38628 ${nCASE}"
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )








