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
#  F: urDEL-signifcant genes, not constrained, near pheno-matched GWAS hit

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










