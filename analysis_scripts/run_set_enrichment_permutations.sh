#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run benchmarking permutation tests
#See set_enrichment_benchmarking.sh for description

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
VF=$1
CNV=$2
pheno=$3
n=$4
size=$5
sd=$6

#####run
for i in $( seq -w 0001 1000 ); do
	echo -e "STARTING TEST ${i}"
	${WRKDIR}/bin/rCNVmap/bin/annoSet_permutation_test.sh -N 10000 \
	-x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
	-o ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt \
	${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
	${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
	${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed.gz \
	${WRKDIR}/data/misc/GRCh37_autosomes.genome
done