#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate organ-level master noncoding annotation tracks

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

#####Run
while read anno; do
  echo ${anno}
  while read file; do
    cat ${file}
  done < <( awk -v tissue=${tissue} -v anno=${anno} \
  '{ if ($1==tissue && $2==anno) print $3 }' \
  ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list ) | \
  sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/${tissue}_MASTER.${anno}.elements.bed
done < <( awk -v tissue=${tissue} '{ if ($1==tissue) print $2 }' \
  ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
  sort | uniq )