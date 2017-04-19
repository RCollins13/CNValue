#!/usr/bin/env bash

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
chr=$1

#####Run
awk -v OFS="\t" '{ if ($4>=1) print $1, $2, $3 }' \
${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bg | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 100 -i - | \
awk -v OFS="\t" '{ if ($3-$2>=200) print $1, $2, $3 }' | sed 's/^chr//g' > \
${WRKDIR}/data/misc/PhyloP/evolutionarilyConserved_PhyloP.chr${chr}.elements.bed