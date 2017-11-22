#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to intersect DD-associated gene lists to derive a high-confidence set

#Master annotation curation code

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Prep directories
if ! [ -e ${WRKDIR}/analysis/DD_geneLists ]; then
	mkdir ${WRKDIR}/analysis/DD_geneLists
fi

###################################
#####Define intermediate gene lists
###################################

#####Statistical genetic association with DDs
#DDD
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list | \
sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/DDD.genes.list
#extTADA (ASD/ID/DD)
for pheno in ASD DD ID; do
  sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/extTADA_${pheno}.genes.list
done | sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/extTADA.genes.list
#Coe 2017 (LGD)
for set in Coe2017_ASD_ID_CHModel_LGD Coe2017_ASD_ID_denovolyzeR_LGD; do
  sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/${set}.genes.list
done | sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/Coe_LGD.genes.list
#Coe 2017 (Missense)
for set in Coe2017_ASD_ID_CHModel_MIS Coe2017_ASD_ID_CHModel_MIS30 Coe2017_ASD_ID_CLUMP_MIS Coe2017_ASD_ID_denovolyzeR_MIS; do
  sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/${set}.genes.list
done | sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/Coe_MIS.genes.list

#####Manually curated DD gene sets
#DECIPHER LOF, GOF, and unknown
for mode in LOF GOF Unknown; do
  cat ${WRKDIR}/data/master_annotations/genelists/DDG2P_*_${mode}.genes.list | \
  sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/DDG2P_${mode}.genes.list
done
#SFARIGene

ASD_SFARIGene.genes.list
ClinGen_haploinsufficient_low_confidence.genes.list
#OMIM
Neurodevelopmental_disorder_HPO_associated
Autism_spectrum_disorder_HPO_associated
Developmental_delay_HPO_associated
Intellectual_disability_HPO_associated

#####Functional & predicted evidence
#GTEx brain high expressors
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/BRAIN_MASTER_INTERSECTION.Highly_Expressed.genes.list | \
sort | uniq | fgrep -wf - \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/GTEx_brainHighExpressors.genes.list
#GO brain, embryonic, and post-embryonic development
for set in GO_brain_development GO_in_utero_embryonic_development GO_post_embryonic_development; do
  sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/${set}.genes.list | \
  sort | uniq | sed 's/\-/_/g' | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/${set}.genes.list
done
ASD_BrainExpressionNetworks.genes.list
ASD_DAWN.genes.list

#####Mutational constraint
#ExAC LoF constrained
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
sort | uniq | fgrep -wf - \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/ExAC_LoF_constrained.genes.list
#ExAC missense constrained
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_missense_constrained.genes.list | \
sort | uniq | fgrep -wf - \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/ExAC_missense_constrained.genes.list
#RVIS intolerant
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/RVIS_intolerant.genes.list | \
sort | uniq | fgrep -wf - \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
grep -e '^[0-9]' | cut -f4 | sort | uniq > ${WRKDIR}/analysis/DD_geneLists/RVIS_intolerant.genes.list

#####









