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
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

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
	${WRKDIR}/bin/rCNVmap/bin/annoSet_permutation_test.sh -N 1000 \
	-x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
	-o ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt \
	${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
	${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
	${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed.gz \
	${WRKDIR}/data/misc/GRCh37_autosomes.genome
done