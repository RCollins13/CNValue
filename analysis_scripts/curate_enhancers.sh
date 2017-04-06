#!/usr/bin/env bash

#TEST

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate enhancers & target genes from EnhancerAtlas

#####Set parameters
WRKDIR=/data/talkowski/Samples/rCNVmap
cd ${WRKDIR}
h37=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.fa
TMPDIR=/scratch/miket/rlc47temp/tmp.files/
module load bedtools/2.22.1
module load samtools/1.3
module load hdf/1.8.14
module load anaconda/4.0.5
module rm gcc-4.4
module load gcc/4.9.0
module load bcftools/1.3.1
PICARD=/data/talkowski/tools/bin/picard-tools-1.137/picard.jar
SFARI_ANNO=/data/talkowski/Samples/SFARI/ASC_analysis/annotations

#####Read arguments
tissue=$1

#####run
for i in $( seq -w 0001 0100 | paste -s ); do
  if [ ${W} == 0 ]; then
  	if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_exons_n${n}_i${i}.txt ]; then 
  	  fgrep -v "#" ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_exons_n${n}_i${i}.txt | awk '{ print $NF }'
  	fi
  else
    if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_wholegenes_n${n}_i${i}.txt ]; then 
      fgrep -v "#" ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_wholegenes_n${n}_i${i}.txt | awk '{ print $NF }'
    fi
  fi
done > ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}/${pheno}/${VF}_results_n${n}_W${W}.txt