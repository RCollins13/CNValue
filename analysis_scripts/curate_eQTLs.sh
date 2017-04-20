#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate GTEx eQTLs & target genes

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read argument
tissue=$1

#####Run
ntissue=$( echo ${tissue} | sed 's/\-/_/g' )
zcat ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL/${tissue}_Analysis.v6p.signif_snpgene_pairs.txt.gz | \
cut -f1-2 | sed 's/\./\t/g' | cut -f1-2 | sed '1d' | sed 's/_/\t/g' | awk -v OFS="\t" \
'{ print $1, $2, $2+length($4), $6 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -c 4 -o distinct -i - | \
sed -f ${WRKDIR}/data/master_annotations/gencode/ENSG_to_symbols.sed > \
${WRKDIR}/data/master_annotations/noncoding/eQTLs_${ntissue}.elements.bed