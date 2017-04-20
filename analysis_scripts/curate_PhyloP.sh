#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate PhyloP conserved elements by chromosome

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
chr=$1

#####Run
awk -v OFS="\t" '{ if ($4>=1) print $1, $2, $3 }' \
${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bg | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 100 -i - | \
awk -v OFS="\t" '{ if ($3-$2>=200) print $1, $2, $3 }' | sed 's/^chr//g' > \
${WRKDIR}/data/misc/PhyloP/evolutionarilyConserved_PhyloP.chr${chr}.elements.bed