#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to generate 1,000 random genomic intervals
#Size is normally distributed
#See set_enrichment_benchmarking.sh for description

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
size=$1
sd=$2
n=$3

#####$un
for i in $( seq -w 0001 1000 ); do
  ${WRKDIR}/bin/rCNVmap/bin/simulate_intervals.sh -z -s ${size} -d ${sd} -N ${n} \
   -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
   -o ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed \
   ${WRKDIR}/data/misc/GRCh37_autosomes.genome
done