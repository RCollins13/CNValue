#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Helper script to parallelize filtering of master burden file

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
list=`mktemp`
echo -e "${group}.${CNV}.${filt}.obs_p" > ${list}
${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
-o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R -t 0.0000003818461 \
-o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
-o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
rm ${list}
#Make list of loci (±20kb merge distance & min 20kb size)
ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct -d 20000 \
-i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | \
awk '{ if ($3-$2>20000) print $0 }' > \
${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed

