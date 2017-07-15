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

#####Run geneScore model for all combinations
#Create master output directory
if [ -e ${WRKDIR}/analysis/perGene_burden ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden
fi
mkdir ${WRKDIR}/analysis/perGene_burden
#Run model
while read pheno; do
  #Create output directory for results
  if [ -e ${WRKDIR}/analysis/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perGene_burden/${pheno}
  #Get number of subjects in group
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
           ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${context}_geneScoreModel \
          "${WRKDIR}/bin/rCNVmap/bin/run_geneScore_model.R \
          -o ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
          ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          38628 ${nCASE}"
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )

#####Subselect significant genes per combination
#Create master output directory
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes
#Run model
while read pheno; do
  #Create output directory for results
  if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          #Nominal
          awk -v FS="\t" '{ if ($51<=0.05) print $1 }' \
          ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_nominally_sig.genes.list
          #FDR
          awk -v FS="\t" '{ if ($52<=0.05) print $1 }' \
          ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_FDR_sig.genes.list
          #Bonferroni
          awk -v FS="\t" '{ if ($53<=0.05) print $1 }' \
          ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_Bonferroni_sig.genes.list
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )

#####Get summary table of significant gene counts
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do




#####Test code for enrichment tests vs gene sets
#Get NDD-sig E4 del genes
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/NDD/NDD_DEL_E4_exonic.geneScore_stats.txt > \
${TMPDIR}/NDD_DEL_E4_exonic.sig_genes.txt
#SEIZ E4 DEL
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/SEIZ/SEIZ_DEL_E4_exonic.geneScore_stats.txt > \
${TMPDIR}/SEIZ_DEL_E4_exonic.sig_genes.txt
#SCZ E4 DEL
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/SCZ/SCZ_DEL_E4_exonic.geneScore_stats.txt > \
${TMPDIR}/SCZ_DEL_E4_exonic.sig_genes.txt
#SEIZ E4 DUP
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/SEIZ/SEIZ_DEL_E2_wholegene.geneScore_stats.txt > \
${TMPDIR}/SEIZ_DUP_E2_wholegene.sig_genes.txt
#CARD E4 DEL
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/CARD/CARD_DEL_E4_exonic.geneScore_stats.txt > \
${TMPDIR}/CARD_DEL_E4_exonic.sig_genes.txt
#Get CTRL-sig E4 del genes
awk '{ if ($35<=0.05) print $1 }' \
${WRKDIR}/analysis/perGene_burden/DD/DD_DEL_E4_exonic.geneScore_stats.txt > \
${TMPDIR}/DD_DEL_E4_exonic.sig_genes.txt
#Get autosomal protein-coding genes
grep -ve '^X\|^Y\|^M' \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f4 | sort | uniq > ${TMPDIR}/autosomal_protein_coding_genes.list
#Run tests
/data/talkowski/rlc47/code/ScriptToolbox/bidirectionalEnrichment.sh -n 20 \
${TMPDIR}/DD_DEL_E4_exonic.sig_genes.txt \
${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list \
${TMPDIR}/autosomal_protein_coding_genes.list
/data/talkowski/rlc47/code/ScriptToolbox/bidirectionalEnrichment.sh -n 20 \
${TMPDIR}/DD_DEL_E4_exonic.sig_genes.txt \
${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list \
${TMPDIR}/autosomal_protein_coding_genes.list








