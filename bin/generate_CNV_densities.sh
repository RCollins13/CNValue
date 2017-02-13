#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Generate CNV density bedgraphs for a given analysis group

#Set parameters
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

#Read args
group=$1

#Run
for CNV in DEL DUP CNV; do
  echo ${CNV}
  bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz > \
  ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
  gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
  bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz > \
  ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
  gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
  bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz > \
  ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
  gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
done