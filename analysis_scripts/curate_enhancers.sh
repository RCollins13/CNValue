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
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
tissue=$1

#####run
fgrep ">" ${WRKDIR}/data/misc/EnhancerAtlas/${tissue}.fasta | \
sed -e 's/>//g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/_/\t/g' -e 's/^chr//g' | \
awk -v OFS="\t" '{ if ($4=="") $4="."; print }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - | \
sed -f ${WRKDIR}/data/master_annotations/gencode/ENSG_to_symbols.sed > \
${WRKDIR}/data/master_annotations/noncoding/enhancers_${tissue}.elements.bed