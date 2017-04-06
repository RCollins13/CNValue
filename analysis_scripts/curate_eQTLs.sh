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