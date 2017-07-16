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
for pheno in GERM NEURO NDD PSYCH SOMA CNCR; do
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






