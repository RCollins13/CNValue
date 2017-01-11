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
  col=$( zcat ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz | \
  head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${group}.${CNV}.${filt}.obs_p" ) | cut -f1 )
  bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
  -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
  cut -f1-4,${col} | sort -nk5,5 | head -n1 | cut -f1-4
done < <( bedtools intersect -u -wa \
  -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
  -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
  cut -f1-3 | fgrep -v "#" ) > \
${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
