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

#####Gencode
mkdir ${WRKDIR}/data/master_annotations/gencode
#Download & unpack
cd ${WRKDIR}/data/master_annotations/gencode
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf.gz
#Cut list of all gene symbols
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
awk '{ if ($3=="gene") print $0 }' | sed 's/\;/\n/g' | fgrep "gene_name" | tr -d "\"" | \
awk '{ print $2 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/Gencode_v19_all.genes.list
#Cut list gene symbols for select transcript types
for class in antisense lincRNA miRNA protein_coding pseudogene snoRNA snRNA; do
  echo ${class}
  fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
  awk '{ if ($3=="gene") print $0 }' | fgrep "transcript_type \"${class}\"" | \
  sed 's/\;/\n/g' | fgrep "gene_name" | tr -d "\"" | awk '{ print $2 }' | sort | uniq > \
  ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_${class}.genes.list
done

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
#Constrained genes (pLI > 0.9)
cd ${WRKDIR}/data/misc/
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
awk '{ if ($20>0.9 && $20!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list
#Dispensable genes (pLI < 0.1)
awk '{ if ($20<0.1 && $20!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_dispensable.genes.list
#Highly constrained genes (pLI > 0.9 & LOF obs:exp top half among constrained genes )
nconst=$( sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk '{ if ($20>0.9 && $20!="NA" && $16>0) print $0 }' | wc -l )
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($20>0.9 && $20!="NA" && $16>0) print $13/$16, $2 }' | \
sort -nk1,1 | head -n $((${nconst}/2)) | cut -f2 | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_highly_constrained.genes.list
#Extremely constrained genes (pLI > 0.9 & no LoF mutations observed in ExAC)
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($20>0.9 && $20!="NA" && $16>0 && $13==0) print $2 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_extremely_constrained.genes.list
#Intolerant genes (RVIS top 10% for MAF < 0.1% in any ExAC population)
cd ${WRKDIR}/data/misc/
wget http://genic-intolerance.org/data/RVIS_Unpublished_ExAC_May2015.txt
sed '1d' ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | \
awk '{ if ($11<10) print $5 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/RVIS_intolerant.genes.list
#Highly intolerant genes (RVIS top 5% for MAF < 0.1% in any ExAC population)
sed '1d' ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | \
awk '{ if ($11<5) print $5 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/RVIS_highly_intolerant.genes.list
#Extremely intolerant genes (RVIS top 1% for MAF < 0.1% in any ExAC population)
sed '1d' ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | \
awk '{ if ($11<1) print $5 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/RVIS_extremely_intolerant.genes.list
#OMIM disease-associated genes
cd ${WRKDIR}/data/misc/
git clone https://github.com/macarthur-lab/gene_lists.git
mv gene_lists MacArthur_gene_lists
sed '1d' ${WRKDIR}/data/misc/MacArthur_gene_lists/other_data/omim.full.tsv | \
awk -v OFS="\n" '$4 !~ /NA/ { print $1, $3 }' | sed -e 's/|/\n/g' -e 's/,/\n/g' | \
awk '{ if ($1!="NA") print $0 }' | sort | uniq | sed '/^$/d' > \
${WRKDIR}/data/master_annotations/genelists/OMIM_all.genes.list




