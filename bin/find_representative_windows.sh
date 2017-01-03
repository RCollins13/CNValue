#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Parallelization script for finding peak (most significant) windows per
# final locus per disease comparison

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
group=$1
CNV=$2
filt=$3

#####Run
while read chr start end; do
  size=$((${end}-${start}))
  mod=$( expr ${size} % 100000 )
  case ${mod} in
    0)
      paste <( seq ${start} 100000 $((${end}-100000)) ) \
      <( seq $((${start}+100000)) 100000 ${end} ) | \
      awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      bedtools intersect -wa -r -f 1 -b - \
      -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      cut -f1-4
      ;;
    25000)
      paste <( seq ${start} 100000 $((${end}-125000)) ) \
      <( seq $((${start}+100000)) 100000 $((${end}-25000)) ) | \
      awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      bedtools intersect -wa -r -f 1 -b - \
      -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      cut -f1-4
      ;;
    50000)
      paste <( seq $((${start}-25000)) 100000 $((${end}-75000)) ) \
      <( seq $((${start}+75000)) 100000 $((${end}+25000)) ) | \
      awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      bedtools intersect -wa -r -f 1 -b - \
      -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      cut -f1-4
      ;;
    75000)
      paste <( seq ${start} 100000 $((${end}-75000)) ) \
      <( seq $((${start}+100000)) 100000 $((${end}+25000)) ) | \
      awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      bedtools intersect -wa -r -f 1 -b - \
      -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      cut -f1-4
      ;;
  esac
done < <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
  fgrep -v "#" | cut -f1-3 ) | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
