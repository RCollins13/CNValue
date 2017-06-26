#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure 3

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/figure3 ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure3
fi
mkdir ${WRKDIR}/data/plot_data/figure3

#####Get obs/exp constraint scores for genes from gene sets signif in 
#####0/35, 1/35, and 35/35 pheno groups
#Make lists of genes for each category
for criteria in onePheno allPhenos; do
  while read geneset; do
    awk -v OFS="\t" -v geneset=${geneset} '{ if ($1==geneset) print $2 }' \
    ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.list | xargs -I {} cat {}
  done < <( fgrep -wv "All_protein_coding_genes" \
            ${WRKDIR}/bin/rCNVmap/misc/geneSets_signif_${criteria}.list | \
            fgrep -wv "All_Gencode_genes" ) | \
  sort | uniq | \
  fgrep -wf ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list > \
  ${TMPDIR}/${criteria}.genes.list
done
for criteria in onePheno allPhenos; do
  cat ${TMPDIR}/${criteria}.genes.list
done | sort | uniq > ${TMPDIR}/anyPhenos.genes.list
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | \
fgrep -wvf <( cat ${TMPDIR}/onePheno.genes.list \
                  ${TMPDIR}/allPhenos.genes.list | \
              sort | uniq | sed 's/\-/_/g' ) > \
${TMPDIR}/noPhenos.genes.list
#Collect ExAC obs:exp z-scores
paste <( cut -f2 ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
sed '1d' | sed 's/\-/_/g' ) \
<( cut -f19 ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sed '1d' ) > \
${TMPDIR}/ExAC_LOF_z.list
#Extract relevant z-scores per list
for criteria in noPhenos onePheno allPhenos anyPhenos; do
  fgrep -wf ${TMPDIR}/${criteria}.genes.list ${TMPDIR}/ExAC_LOF_z.list | \
  cut -f2 > ${WRKDIR}/data/plot_data/figure3/ExAC_LoF_z.${criteria}.txt
done
cut -f19 ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
sed '1d' > ${WRKDIR}/data/plot_data/figure3/ExAC_LoF_z.allGenes.txt





