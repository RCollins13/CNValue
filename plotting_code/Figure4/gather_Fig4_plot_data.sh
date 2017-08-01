#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure 4

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/figure4 ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure4
fi
mkdir ${WRKDIR}/data/plot_data/figure4

#####Copy geneScore results for GERM/NEURO/NDD/PSYCH/SOMA/CNCR, exonic/wholegene, E3/E4/N1, CNV/DEL/DUP
#Initialize directory
if [ -e ${WRKDIR}/data/plot_data/figure4/geneScore_results ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure4/geneScore_results
fi
mkdir ${WRKDIR}/data/plot_data/figure4/geneScore_results
#Copy all files
for pheno in GERM NEURO NDD ASD PSYCH SOMA HEAD CNCR; do
  for CNV in CNV DEL DUP; do
    for VF in E3 E4 N1; do
      for context in exonic wholegene; do
        cp ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
        ${WRKDIR}/data/plot_data/figure4/geneScore_results/
      done
    done
  done
done

#####Cut ExAC lof obs/exp, Z-scores, and pLI for correlation vs CNV Z-score
#Constraint
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($16>0) print $2, $13/$16, $19, $20 }' > \
${WRKDIR}/data/plot_data/figure4/ExAC_LoF_constraint.txt
#RVIS (0.01%)
sed '1d' ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | \
awk -v OFS="\t" '{ print $1, $(NF-1), $NF }' > \
${WRKDIR}/data/plot_data/figure4/ExAC_RVIS.txt

#####Gather data per phenotype for gene set enrichment plots
#Make universal list of genes tested
fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/GERM/GERM_CNV_E2_exonic.geneScore_stats.txt | \
cut -f1 > ${TMPDIR}/all_tested_genes.list
#Create directories
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/
#Submit data collection - merged master pheno groups
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/geneSet_enrichment_comparison_sets.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/${CNV}_${VF}_${context}_${sig}.comparisons.txt"
      done
    done
  done
done
#All original pheno groups, all gene sets
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons_allPhenos_allGenes -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.all_phenos.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/${CNV}_${VF}_${context}_${sig}.comparisons.txt"
      done
    done
  done
done
#All original pheno groups, all gene sets, exclude constrained genes
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons_allPhenos_allGenes_noConstrained -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.all_phenos.wExclusion.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/${CNV}_${VF}_${context}_${sig}.comparisons.txt \
        ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list"
      done
    done
  done
done
#Copy to plot data directory
if [ -e ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons
fi
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/subset/
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets_noConstrained/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/subset/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets_noConstrained/

#####Copy significant genes to plotting data directory
#Make directory
if [ -e ${WRKDIR}/data/plot_data/signif_genes_unique ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_genes_unique
fi
mkdir ${WRKDIR}/data/plot_data/signif_genes_unique
cp /data/talkowski/Samples/rCNVmap/analysis/perGene_burden/signif_genes/merged/*unique.genes.list \
${WRKDIR}/data/plot_data/signif_genes_unique/

#####Get union set of NDD and CNCR E4 exonic Bonf sig genes
VF="E4"
for CNV in DEL DUP; do
  for pheno in CNCR NDD; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list
  done
done | sort | uniq > \
${TMPDIR}/NDD_CNCR_DEL_DUP_Bonferroni.genes.list
while read gene; do
  for dummy in 1; do
    echo ${gene}
    for pheno in NDD CNCR; do
      for CNV in DEL DUP; do
        fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list | wc -l
      done
    done
  done | paste -s
done < ${TMPDIR}/NDD_CNCR_DEL_DUP_Bonferroni.genes.list > \
${WRKDIR}/data/plot_data/figure4/NDD_CNCR_PPI_CNVmembership.txt
for CNV in CNV DEL DUP; do
  for pheno in CNCR NDD; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list
  done
done | sort | uniq > \
${TMPDIR}/NDD_CNCR_CNV_DEL_DUP_Bonferroni.genes.list
while read gene; do
  for dummy in 1; do
    echo ${gene}
    for pheno in NDD CNCR; do
      for CNV in CNV DEL DUP; do
        fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list | wc -l
      done
    done
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | wc -l
  done | paste -s
done < ${TMPDIR}/NDD_CNCR_CNV_DEL_DUP_Bonferroni.genes.list > \
${WRKDIR}/data/plot_data/figure4/NDD_CNCR_PPI_CNVmembership_CNV_DEL_DUP.txt

#####Get counts of significant genes (DEL/DUP) per phenotype
VF=E4
context=exonic
sig=Bonferroni
for dummy in 1; do
  #Control
  echo "CTRL"
  for CNV in DEL DUP; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
  done
  #Affected pheno groups
  while read pheno; do
    echo "${pheno}"
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
    done
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
            fgrep -v CTRL | cut -f1 )
  #Affected - unions
  for group in GERM NEURO NDD PSYCH SOMA CNCR; do
    echo "ALL_${group}_UNION"
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
    done
  done
done | paste - - -






