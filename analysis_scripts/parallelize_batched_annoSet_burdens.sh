#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to parallelize annotation set CNV burden tests by phenotype-CNV combo

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
pheno=$1
CNV=$2

#####Run
for filt in all coding haplosufficient noncoding intergenic; do
  for VF in E2 E3 E4 N1; do
    echo "\n\n${pheno} ${CNV} ${filt} ${CNV}\n\n"
    #Submit batch job
    ${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
      -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
      -p ${pheno}_${CNV}_${filt}_${VF} \
      -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list \
      ${WRKDIR}/data/misc/GRCh37_autosomes.genome
  done
done