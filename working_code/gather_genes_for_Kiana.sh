#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather list of genes of interest for Kiana's literature search & manual curation

#Gene criteria:
#  A+B: Top cancer-germline pleiotropic genes with opposing DEL/DUP biases
#  C: NDD urDEL-signifcant genes not in DDD 2017 LOF exome list
#  D: urDUP-significant genes, not constrained, highly expressed in one tissue and silent in another tissue
#  E: urDEL- or urDUP-signifcant genes, not constrained, near pheno-matched GWAS hit
#  F: whole-gene urDUP-significant genes, germline, ranked by Z-score

####################
#####Load parameters
####################
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

##########
#####Setup
##########
if [ -e ${TMPDIR}/genes_for_KM ]; then
  rm -rf ${TMPDIR}/genes_for_KM
fi
mkdir ${TMPDIR}/genes_for_KM

#####################################
#####Get lists of genes: Groups A & B
#####################################
#Group A: cancer-germ pleio genes with DEL in GERM & DUP in CNCR
#Group B: cancer-germ pleio genes with DUP in GERM & DEL in CNCR
#Criteria: must be in >20% of GERM & CNCR pheno groups for DEL-DUP or DUP-DEL
VF=E4
context=exonic
sig=Bonferroni
thresh=0.20
#Pleiotropic GERM & CNCR genes
for group in GERM CNCR; do
  nGroups=$( awk -v group=${group} '{ if ($2==group) print $1 }' \
              ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | wc -l )
  for CNV in DEL DUP; do
    while read pheno; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
              ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
    sort | uniq -c | awk -v OFS="\t" -v thresh=${thresh} -v nGroups=${nGroups} \
    '{ if ($1/nGroups>=thresh) print $2 }' > ${TMPDIR}/${group}_${CNV}_pleio_${thresh}pctMin.genes.list
  done
done
#Cross-compare
fgrep -wf ${TMPDIR}/GERM_DEL_pleio_${thresh}pctMin.genes.list \
${TMPDIR}/CNCR_DUP_pleio_${thresh}pctMin.genes.list > \
${TMPDIR}/genes_for_KM/GroupA.genes.list
fgrep -wf ${TMPDIR}/GERM_DUP_pleio_${thresh}pctMin.genes.list \
${TMPDIR}/CNCR_DEL_pleio_${thresh}pctMin.genes.list > \
${TMPDIR}/genes_for_KM/GroupB.genes.list

################################
#####Get lists of genes: Group C
################################
#Group C: constrained, DEL-sig NDD genes not in DDD 2017 list
VF=E4
context=exonic
sig=Bonferroni
fgrep -wvf ${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list \
${WRKDIR}/analysis/perGene_burden/signif_genes/merged/NDD_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | \
fgrep -wf ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list > \
${TMPDIR}/genes_for_KM/GroupC.genes.list

################################
#####Get lists of genes: Group D
################################
#Group D: not constrained, DUP-sig GERM genes, highly expressed in one tissue and not expressed in another
VF=E4
context=exonic
sig=Bonferroni
#Get list of unconstrained high expressors in any tissue
cat ${WRKDIR}/data/master_annotations/genelists/*highly_expressed* | sort | uniq | \
fgrep -wvf ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
sed 's/\-/_/g' | fgrep -wf \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list ) > \
${TMPDIR}/highlyExpressed_notConstrained.genes.list
#Get list of unconstrained silent genes in any tissue
cat ${WRKDIR}/data/master_annotations/genelists/*not_expressed* | sort | uniq | \
fgrep -wvf ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
sed 's/\-/_/g' | fgrep -wf \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list ) > \
${TMPDIR}/notExpressed_notConstrained.genes.list
#Get union of previous two lists
fgrep -wf ${TMPDIR}/highlyExpressed_notConstrained.genes.list \
${TMPDIR}/notExpressed_notConstrained.genes.list > \
${TMPDIR}/notExpressedAndHighlyExpressed_notConstrained.genes.list
#Extract group D genes
group=GERM
CNV=DUP
thresh=0.20
while read pheno; do
  cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
          ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
sed 's/\-/_/g' | sort | uniq -c | awk -v thresh=${thresh} \
'{ if ($1/22>=thresh) print $2 }' | fgrep -wf \
${TMPDIR}/notExpressedAndHighlyExpressed_notConstrained.genes.list > \
${TMPDIR}/genes_for_KM/GroupD.genes.list

################################
#####Get lists of genes: Group E
################################
#Group E: non-constrained, DEL or DUP significant genes, near (±250kb) GWAS hits from relevant phenotype
VF=E4
context=exonic
sig=Bonferroni
dist=250000
#Get list
while read pheno; do
  for CNV in DEL DUP; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | \
    sed 's/\-/_/g' | fgrep -wvf \
    <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | \
    fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
    awk -v OFS="\t" -v dist=${dist} '{ print $1, $2-dist, $3+dist, $4 }' | \
    awk -v OFS="\t" '{ if ($2<1) $2=1; print }' | bedtools intersect -wb -b - \
    -a ${WRKDIR}/data/master_annotations/noncoding/GWAS_loci_NDD.elements.bed
    done | awk '{ print $NF }' | sort | uniq
done < <( sed '1,2d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | sort | uniq > \
${TMPDIR}/genes_for_KM/GroupE.genes.list

################################
#####Get lists of genes: Group F
################################
#Group F: whole-gene urDUP-significant genes, germline, ranked by Z-score
CNV=DUP
VF=E4
sig=Bonferroni
context=wholegene
while read pheno; do
    fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
    awk -v OFS="\t" -v pheno=${pheno} '{ if ($41>=3) print $1, pheno, $41, $26 }'
done < <( awk '{ if ($2=="GERM") print $1 }' \
          ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort -nrk3,3 | cut -f1 | head -n200 | awk '!x[$0]++' | head -n20 > \
  ${TMPDIR}/genes_for_KM/GroupF.genes.list

###########################################
#####Gather metadata for all genes per list
###########################################
VF=E4
sig=Bonferroni
context=exonic
#Write header to metadata file
paste <( echo -e "gene\tcriteria" ) \
<( sed '1,2d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 | paste -s ) \
<( echo -e "constrained\tgeneSets" ) > \
${WRKDIR}/data/plot_data/topCandidateGenes_forKM.annotated.txt
#Collect metadata and write to file (groups A-E)
while read gene; do
  for dummy in 1; do
    #Print gene symbol
    echo ${gene}
    #Get criteria list
    for set in A B C D E F; do
      if [ $(fgrep -w ${gene} ${TMPDIR}/genes_for_KM/Group${set}.genes.list | wc -l) -gt 0 ]; then
        echo ${set}
      fi
    done | paste -s -d,
    #Iterate over phenotypes and report significant associations per group
    while read pheno; do
      DEL=$( fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      DUP=$( fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      if [ ${DEL} -gt 0 ]; then
        if [ ${DUP} -gt 0 ]; then
          echo "BOTH"
        else
          echo "DEL"
        fi
      else
        if [ ${DUP} -gt 0 ]; then
          echo "DUP"
        else
          echo "NOT"
        fi
      fi
    done < <( sed '1,2d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 )
    #Constraint
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | wc -l
    #Other gene sets
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/* | \
    sed -e 's/\//\t/g' -e 's/\:/\t/g' | awk '{ print $(NF-1) }' | \
    sed 's/\.genes\.list//g' | paste -s -d,
  done | paste -s
done < <( cat ${TMPDIR}/genes_for_KM/GroupA.genes.list \
              ${TMPDIR}/genes_for_KM/GroupB.genes.list \
              ${TMPDIR}/genes_for_KM/GroupC.genes.list \
              ${TMPDIR}/genes_for_KM/GroupD.genes.list \
              ${TMPDIR}/genes_for_KM/GroupE.genes.list | sort | uniq ) >> \
${WRKDIR}/data/plot_data/topCandidateGenes_forKM.annotated.txt
#Collect metadata and write to file (group F)
context=wholegene
while read gene; do
  for dummy in 1; do
    #Print gene symbol
    echo ${gene}
    #Get criteria list
    for set in A B C D E F; do
      if [ $(fgrep -w ${gene} ${TMPDIR}/genes_for_KM/Group${set}.genes.list | wc -l) -gt 0 ]; then
        echo ${set}
      fi
    done | paste -s -d,
    #Iterate over phenotypes and report significant associations per group
    while read pheno; do
      DEL=$( fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      DUP=$( fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      if [ ${DEL} -gt 0 ]; then
        if [ ${DUP} -gt 0 ]; then
          echo "BOTH"
        else
          echo "DEL"
        fi
      else
        if [ ${DUP} -gt 0 ]; then
          echo "DUP"
        else
          echo "NOT"
        fi
      fi
    done < <( sed '1,2d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 )
    #Constraint
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | wc -l
    #Other gene sets
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/* | \
    sed -e 's/\//\t/g' -e 's/\:/\t/g' | awk '{ print $(NF-1) }' | \
    sed 's/\.genes\.list//g' | paste -s -d,
  done | paste -s
done < ${TMPDIR}/genes_for_KM/GroupF.genes.list >> \
${WRKDIR}/data/plot_data/topCandidateGenes_forKM.annotated.txt








