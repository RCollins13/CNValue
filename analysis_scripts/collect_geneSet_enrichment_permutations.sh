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
W=$5

#####run
for i in $( seq -w 0001 1000 | paste -s ); do
  if [ ${W} == 0 ]; then
  	if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_exons_n${n}_i${i}.txt ]; then 
  	  fgrep -v "#" ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_exons_n${n}_i${i}.txt | awk '{ print $NF }'
  	fi
  else
    if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_wholegenes_n${n}_i${i}.txt ]; then 
      fgrep -v "#" ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}/results_wholegenes_n${n}_i${i}.txt | awk '{ print $NF }'
    fi
  fi
done > ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}/${pheno}/${VF}_results_n${n}_W${W}.txt