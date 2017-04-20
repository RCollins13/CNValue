#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to collect results from benchmarking permutation tests

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
for i in $( seq -w 0001 1000 | paste -s ); do 
	if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt ]; then 
	  fgrep -v "#" ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt | awk '{ print $NF }'
	fi
done > ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}/${pheno}/${VF}_results_${size}bp_${sd}bp_x${n}.txt