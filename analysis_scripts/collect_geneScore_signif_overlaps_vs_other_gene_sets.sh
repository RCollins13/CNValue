#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects gene count overlap of autosomal protein-coding genes from unique signif
# per-gene tests (geneScore) vs a supplied set of comparison gene lists

# Automatically runs for CTRL, GERM, NEURO, NDD, PSYCH, SOMA, CNCR (in order)

#Read arguments
CNV=$1
VF=$2
context=$3
signif=$4
complist=$5 #two columns: gene list name & full path. Tab-delimmed
universe=$6 #list of genes considered in test
outfile=$7

#Write header
echo -e "#geneset\tall_genes_tested\tbaseline_count\tall_CTRL\tCTRL_count\t\
all_GERM\tGERM_count\tall_NEURO\tNEURO_count\tall_NDD\tNDD_count\t\
all_PSYCH\tPSYCH_count\tall_SOMA\tSOMA_count\tall_CNCR\tCNCR_count" > ${outfile}

#Iterate over all comparison lists
while read ID list; do
  for dummy in 1; do
    echo ${ID}
    #Get baseline number of all genes
    fgrep -v "#" ${universe} | wc -l
    #Get number of genes in list in universe
    fgrep -wf <( sed 's/\-/_/g' ${list} ) <( sed 's/\-/_/g' ${universe} ) | wc -l
    #Iterate over all phenotypes
    for pheno in CTRL GERM NEURO NDD PSYCH SOMA CNCR; do
      #Get number of genes significant in phenotype
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${pheno}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
      #Get number of genes overlapping comparison list
      fgrep -wf \
      <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${pheno}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list ) \
      <( sed 's/\-/_/g' ${list} ) | wc -l
    done 
  done | paste -s
done < ${complist} >> ${outfile}