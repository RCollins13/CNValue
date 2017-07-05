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

#####Collect constraint decile burden test results
#Prepare directory
mkdir ${WRKDIR}/data/plot_data/figure3/constraint_deciles
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for filt in all; do
      for context in exonic; do
        for collection in effectSize pValue lowerCI upperCI zScore; do
          bsub -q short -sla miket_sc -u nobody \
          -J ${CNV}_${VF}_${filt}_${context}_${collection} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneSet_burdens.constraint_deciles.sh \
          ${CNV} ${VF} ${filt} ${context} ${collection}"
        done
      done
    done
  done
done

#####Collect constraint quintile burden test results
#Prepare directory
mkdir ${WRKDIR}/data/plot_data/figure3/constraint_quintiles
#Grouped by CNV type, VF, and CNV filter. MxN matrix; M: phenos, N: annotation sets
#One matrix of p-values, one matrix of effect sizes, and two matrices of 
# confidence interval bounds per set of filters
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for filt in all; do
      for context in exonic; do
        for collection in effectSize pValue lowerCI upperCI zScore; do
          bsub -q short -sla miket_sc -u nobody \
          -J ${CNV}_${VF}_${filt}_${context}_${collection} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneSet_burdens.constraint_quintiles.sh \
          ${CNV} ${VF} ${filt} ${context} ${collection}"
        done
      done
    done
  done
done

#####Compute similarity between gene sets
#Write header
paste <( echo "SET" ) \
      <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list | paste -s ) > \
${WRKDIR}/data/plot_data/figure3/gene_set_overlaps.matrix.txt
#Iterate over all gene sets
while read nameA setA; do
  for dummy in 1; do
    echo ${nameA}
    n_setA=$( cat ${setA} | wc -l )
    while read nameB setB; do
      if [ -z ${setB} ]; then
        sed 's/\-/_/g' ${setB} | fgrep -wf - \
        <( sed 's/\-/_/g' ${setA} ) | wc -l | \
        awk -v OFS="\t" -v n_setA=${n_setA} '{ print $1/n_setA }'
      else
        echo "NA"
      fi
    done < ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list
  done | paste -s
done < ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list >> \
${WRKDIR}/data/plot_data/figure3/gene_set_overlaps.matrix.txt












