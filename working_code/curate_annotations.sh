#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master annotation curation code

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

#####Prepare annotation directory tree
mkdir ${WRKDIR}/data/master_annotations
mkdir ${WRKDIR}/data/master_annotations/genelists

#####GENE LISTS
#GTEx tissue-specific expressed genes
cp ${SFARI_ANNO}/genelists/GTEx*specific.genes.list \
${WRKDIR}/data/master_annotations/genelists/
while read tissue; do
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list | \
  sort | uniq > ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list2
  mv ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list2 \
  ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list
done < <( l ${WRKDIR}/data/master_annotations/genelists/GTEx_*_specific.genes.list | \
 sed -e 's/GTEx_/\t/g' -e 's/_specific/\t/g' | awk '{ print $(NF-1) }' )
while read tissue; do
  echo "$( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" )-specifically expressed genes"
  echo "Genes preferentially expressed in $( echo ${tissue} | sed 's/_/\ /g' ) from GTEx"
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list | wc -l
  readlink -f ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_specific.genes.list
done < <( l ${WRKDIR}/data/master_annotations/genelists/GTEx_*_specific.genes.list | \
 sed -e 's/GTEx_/\t/g' -e 's/_specific/\t/g' | awk '{ print $(NF-1) }' ) | paste - - - -
#GTEx highly expressed genes (99th percentile)
while read tissue; do
  echo ${tissue}
  Rscript -e "x <- read.table(\"${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz\",\
    header=F); high <- unique(as.character(x[which(x[,5]>=quantile(x[,5],0.99)),4]));\
    write.table(high,\"${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_highly_expressed.genes.list\",\
      col.names=F,row.names=F,quote=F)"
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 )
while read tissue; do
  echo "Highly $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" )-expressed genes"
  echo "Genes highly expressed in $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" ) from GTEx"
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_highly_expressed.genes.list | wc -l
  readlink -f ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_highly_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 ) | paste - - - -
#GTEx very highly expressed genes (99th percentile)
while read tissue; do
  echo ${tissue}
  Rscript -e "x <- read.table(\"${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz\",\
    header=F); high <- unique(as.character(x[which(x[,5]>=quantile(x[,5],0.999)),4]));\
    write.table(high,\"${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_very_highly_expressed.genes.list\",\
      col.names=F,row.names=F,quote=F)"
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 )
while read tissue; do
  echo "Very highly $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" )-expressed genes"
  echo "Genes very highly expressed in $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" ) from GTEx"
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_very_highly_expressed.genes.list | wc -l
  readlink -f ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_very_highly_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 ) | paste - - - -
#GTEx not expressed genes (RPKM < 1)
while read tissue; do
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | \
  awk '{ if ($5<1) print $4 }' | sort | uniq > \
  ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_not_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 )
while read tissue; do
  echo "Not $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" )-expressed genes"
  echo "Genes not expressed in $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" ) from GTEx"
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_not_expressed.genes.list | wc -l
  readlink -f ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_not_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 ) | paste - - - -















