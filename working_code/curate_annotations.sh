#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master annotation curation code
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Prepare annotation directory tree
mkdir ${WRKDIR}/data/master_annotations
mkdir ${WRKDIR}/data/master_annotations/genelists

##########################
##########################
#         GENCODE        #
##########################
##########################
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
#Create master bed file of all exons
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" '{ if ($3=="exon") print $1, $4, $5, $10 }' | \
sed 's/\;/\t/g' | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.all.bed
#Create master bed file of protein-coding exons
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.all.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed
#Create master bed file of protein-coding exons from not-clearly-haplosufficient genes
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_notHaplosufficient.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed
#Get number of elements, median element size, and full path for all gencode exon tracks
for filter in all protein_coding notHaplosufficient; do
  #Count of elements
  cat ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed | wc -l
  #Median element size
  awk '{ print $3-$2 }' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed | \
  sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Count of elements (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed | wc -l
  #Median element size (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed | \
  awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed
done | paste - - - - -
#Create master bed file of gene boundaries
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" '{ if ($3=="gene") print $1, $4, $5, $10 }' | \
sed 's/\;/\t/g' | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed
#Create master bed file of gene boundaries from protein-coding genes
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed
#Create master bed file of gene boundaries from protein-coding, not-clearly-haplosufficient genes
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_notHaplosufficient.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.notHaplosufficient.bed
#Get number of elements, median element size, and full path for all gencode gene tracks
for filter in all protein_coding notHaplosufficient; do
  #Count of elements
  cat ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed | wc -l
  #Median element size
  awk '{ print $3-$2 }' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed | \
  sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Count of elements (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed | wc -l
  #Median element size (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed | \
  awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed
done | paste - - - - -
#Create master bed file of gene promoters
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" '{ if ($3=="gene") print $1, $4, $5, $7, $10 }' | \
sed 's/\;/\t/g' | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4, $5 }' | tr -d "\"" | \
sed 's/^chr//g' | awk -v OFS="\t" '{ if ($4=="+") print $1, $2-5000, $2, $4, $5; else print $1, $3, $3+5000, $4, $5 }' | \
awk -v OFS="\t" '{ if ($2<1) $2=1; print }' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.all.bed
#Create master bed file of gene promoters from protein-coding genes
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.all.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.protein_coding.bed
#Create master bed file of gene boundaries from protein-coding, not-clearly-haplosufficient genes
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_notHaplosufficient.genes.list | \
fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.protein_coding.bed ) | \
sed 's/_/\-/g' > ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.notHaplosufficient.bed
#Get number of elements, median element size, and full path for all gencode promoters
for filter in all protein_coding notHaplosufficient; do
  #Count of elements
  cat ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | wc -l
  #Median element size
  awk '{ print $3-$2 }' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | \
  sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Count of elements (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | wc -l
  #Median element size (autosomal)
  grep -e '^[1-9]' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | \
  awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed
done | paste - - - - -
#Get number of elements, median element size, and full path for all gencode promoter tracks
for filter in all protein_coding notHaplosufficient; do
  #Count of elements
  cat ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | wc -l
  #Median element size
  awk '{ print $3-$2 }' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed | \
  sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.${filter}.bed
done | paste - - -
#Create master sed file for transforming ENSG* IDs into gene symbols
awk -v FS="\t" '{ if ($3=="gene") print $0 }' \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed -e 's/gene_id/\t/g' -e 's/gene_name/\t/g' -e 's/\;\ transcript_id/\t/g' -e 's/\;\ transcript_type/\t/g' | \
awk -v OFS="\t" -v FS="\t" '{ print $10, $12 }' | tr -d "\"" | sed 's/^ //g' | \
sort | uniq > ${TMPDIR}/ENSG_to_symbols.tmp
paste <( cut -f1 ${TMPDIR}/ENSG_to_symbols.tmp | cut -f1 -d\. ) \
<( cut -f2 ${TMPDIR}/ENSG_to_symbols.tmp | sed 's/\-/\\\-/g' ) | \
awk '{ print "s/"$1"/"$2"/g" }' > \
${WRKDIR}/data/master_annotations/gencode/ENSG_to_symbols.sed




##########################
##########################
#       GENE LISTS       #
##########################
##########################
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
#GTEx lowly expressed genes (RPKM < 1)
while read tissue; do
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | \
  awk '{ if ($5<1) print $4 }' | sort | uniq > \
  ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_lowly_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 )
while read tissue; do
  echo "Lowly $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" )-expressed genes"
  echo "Genes lowly expressed in $( echo ${tissue} | sed 's/_/\ /g' | tr "A-Z" "a-z" ) from GTEx"
  cat ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_lowly_expressed.genes.list | wc -l
  readlink -f ${WRKDIR}/data/master_annotations/genelists/GTEx_${tissue}_lowly_expressed.genes.list
done < <( l ${SFARI_ANNO}/misc/GTEx_expression/*bed.gz | fgrep -v highExpressor | \
sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\./\t/g' | cut -f2 ) | paste - - - -
#GTEx not expressed genes (median RPKM=0)
while read tissue; do
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | \
  awk '{ if ($5==0) print $4 }' | sort | uniq > \
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
#Haplosufficient genes (pLI < 0.1)
awk '{ if ($20<0.1 && $20!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_haplosufficient.genes.list
#Missense constrained genes (misZ >= 3.09)
awk '{ if ($18>=3.09 && $18!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sed '1d' | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/ExAC_missense_constrained.genes.list
#Not haplosufficient genes (pLI > 0.1)
awk '{ if ($20>=0.1 && $20!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_notHaplosufficient.genes.list
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
#Constraint deciles (based on obs:exp tranches; n=15,916 genes in constraint file after excluding those with obs:exp=0)
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($13>0) print $2, $13/$16 }' | sort -nk2,2 | cut -f1 | \
split -l 1592 -d -a 1 /dev/stdin ${TMPDIR}/constraint_deciles
for i in $( seq 0 9 ); do
  mv ${TMPDIR}/constraint_deciles${i} \
  ${WRKDIR}/data/master_annotations/genelists/ExAC_constraint_decile_$( echo ${i}+1 | bc ).genes.list
done
#Constraint quintiles (based on obs:exp tranches; n=15,916 genes in constraint file after excluding those with obs:exp=0)
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($13>0) print $2, $13/$16 }' | sort -nk2,2 | cut -f1 | \
split -l 3183 -d -a 1 /dev/stdin ${TMPDIR}/constraint_quintiles
for i in $( seq 0 4 ); do
  mv ${TMPDIR}/constraint_quintiles${i} \
  ${WRKDIR}/data/master_annotations/genelists/ExAC_constraint_quintile_$( echo ${i}+1 | bc ).genes.list
done
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
#FDA drug targets
cd ${WRKDIR}/data/misc/
git clone https://github.com/macarthur-lab/gene_lists.git
mv gene_lists MacArthur_gene_lists
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/fda_approved_drug_targets.tsv | \
sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/FDA_drug_targets.genes.list
#Autosomal dominant disease genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/berg_ad.tsv \
${WRKDIR}/data/misc/MacArthur_gene_lists/lists/blekhman_ad.tsv | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/Autosomal_dominant_disease.genes.list
#Autosomal recessive disease genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/berg_ar.tsv \
${WRKDIR}/data/misc/MacArthur_gene_lists/lists/blekhman_ar.tsv | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/Autosomal_recessive_disease.genes.list
#Culture-essential genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/core_essentials_hart.tsv | \
sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/Culture_essential.genes.list
#DNA repair genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/DRG_KangJ.tsv \
${WRKDIR}/data/misc/MacArthur_gene_lists/lists/DRG_WoodRD.tsv | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/DNA_repair.genes.list
#Phenix + UberPheno gene-to-phenotype linker
#Note: requires Claire's curation & merging of UberPheno & Phenix
#Note: also requires Kiana's parsing of Claire's file into per-phenotype lists
while read upper; do
  lower=$( echo ${upper} | tr "A-Z" "a-z" )
  cut -f2 ${WRKDIR}/data/HPO_map/Genes-HPO-DiseaseStates/HPO_Genes/${lower} | \
  sort | uniq > ${WRKDIR}/data/master_annotations/genelists/${upper}_HPO_associated.genes.list
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
  fgrep -v CTRL | fgrep -v UNK | cut -f4 )
#GWAS nearest genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/gwascatalog.tsv | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/GWAS_nearest.genes.list
#ClinGen haploinsufficient genes
cd ${WRKDIR}/data/misc
wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_gene_curation_list.tsv
fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
awk -v FS="\t" '{ if ($5==3) print $1 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_high_confidence.genes.list
fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
awk -v FS="\t" '{ if ($5==3 || $5==2) print $1 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_medium_confidence.genes.list
fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
awk -v FS="\t" '{ if ($5==3 || $5==2 || $5==1) print $1 }' | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_low_confidence.genes.list
# fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
# awk -v FS="\t" '{ if ($10==3) print $1 }' | sort | uniq > \
# ${WRKDIR}/data/master_annotations/genelists/ClinGen_triploinsufficient_high_confidence.genes.list
# fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
# awk -v FS="\t" '{ if ($10==3 || $10==2) print $1 }' | sort | uniq > \
# ${WRKDIR}/data/master_annotations/genelists/ClinGen_triploinsufficient_medium_confidence.genes.list
# fgrep -v "#" ${WRKDIR}/data/misc/ClinGen_gene_curation_list.tsv | \
# awk -v FS="\t" '{ if ($10==3 || $10==2 || $10==1) print $1 }' | sort | uniq > \
# ${WRKDIR}/data/master_annotations/genelists/ClinGen_triploinsufficient_low_confidence.genes.list
#Olfactory receptors
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/olfactory_receptors.tsv | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/Olfactory_receptors.genes.list
#ClinVar disease-associated genes
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/clinvar_path_likelypath.tsv | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/ClinVar_disease_associated.genes.list
#Kinases
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/kinases.tsv | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/Kinases.genes.list
#G-protein coupled receptors
cat ${WRKDIR}/data/misc/MacArthur_gene_lists/lists/gpcr.tsv | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/GPCRs.genes.list
#UW cancer & immunodeficiency panels -- manually curated, see:
# http://depts.washington.edu/labweb/Divisions/MolDiag/MolDiagGen/index.htm
#Gene lists split by GO annotation
mkdir ${WRKDIR}/data/misc/GO
cd ${WRKDIR}/data/misc/GO
wget http://geneontology.org/gene-associations/goa_human.gaf.gz
gunzip ${WRKDIR}/data/misc/GO/goa_human.gaf.gz
wget http://purl.obolibrary.org/obo/go/go-basic.obo
grep -e '^id\:\ GO\:' -e '^name\:\ ' ${WRKDIR}/data/misc/GO/go-basic.obo | \
paste - - | grep -e '^id' | sed 's/^id\:\ //g' | sed 's/name\:\ //g' | \
sed -e 's/\-/_/g' -e 's/\./_/g' -e 's/\ /_/g' -e 's/\,/_/g' > \
${WRKDIR}/data/misc/GO/all_go_terms.txt
mkdir ${WRKDIR}/data/misc/GO/all_GO_terms_splits
split -a 4 -d -l 100 ${WRKDIR}/data/misc/GO/all_go_terms.txt \
${WRKDIR}/data/misc/GO/all_GO_terms_splits/GO_term_split.
mkdir ${WRKDIR}/data/misc/GO/all_GO_gene_lists
for i in $( seq -w 0000 0465 ); do
  bsub -q short -sla miket_sc -J GO_curation_${i} -u nobody \
  "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_GO_genes.sh \
   ${WRKDIR}/data/misc/GO/all_GO_terms_splits/GO_term_split.${i}"
done
wc -l ${WRKDIR}/data/misc/GO/all_GO_gene_lists/*list | awk -v OFS="\t" \
'{ if ($1>=50 && $1<=5000) print $1, $2 }' | sort -nrk1,1 > \
${WRKDIR}/data/misc/GO/eligible_GO_term_gene_lists.txt
while read GO descrip group; do
  nocolon=$( echo ${GO} | sed 's/\:/_/g' )
  ldescrip=$( echo ${descrip} | tr "A-Z" "a-z" )
  cat ${WRKDIR}/data/misc/GO/all_GO_gene_lists/${nocolon}_${descrip}.genes.list | \
  sort | uniq > ${WRKDIR}/data/master_annotations/genelists/GO_${ldescrip}.genes.list
done < ${WRKDIR}/bin/rCNVmap/misc/GO_term_gene_sets.list
while read group; do
  while read GO descrip skip; do
    nocolon=$( echo ${GO} | sed 's/\:/_/g' )
    cat ${WRKDIR}/data/misc/GO/all_GO_gene_lists/${nocolon}_${descrip}.genes.list
  done < <( awk -v group=${group} '{ if ($3==group) print $0 }' \
    ${WRKDIR}/bin/rCNVmap/misc/GO_term_gene_sets.list ) | sort | uniq > \
    ${WRKDIR}/data/master_annotations/genelists/GO_${group}_all.genes.list
done < <( cut -f3 ${WRKDIR}/bin/rCNVmap/misc/GO_term_gene_sets.list | \
sort | uniq | fgrep -v Other )
#DDD genes
cat ${SFARI_ANNO}/genelists/DDD_2016.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list
#CHD8 targets
cat ${SFARI_ANNO}/genelists/Sugathan2014_CHD8.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/CHD8_targets.genes.list
#TADA genes
for cutoff in 0.1 0.3; do
  cat ${SFARI_ANNO}/genelists/Sanders2015_TADAq_${cutoff}.genes.list | sort | uniq > \
  ${WRKDIR}/data/master_annotations/genelists/ASD_TADA_q${cutoff}.genes.list
done
#DAWN genes
cat ${SFARI_ANNO}/genelists/Liu2014_DAWN.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ASD_DAWN.genes.list
#Turner exome DNM ASD genes 
cat ${SFARI_ANNO}/genelists/ASD_Turner57.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ASD_WES_DNM.genes.list
#Note: manually curated SFARIGene & Brain Expression Network ASD lists
#RBFOX1 targets
cat ${SFARI_ANNO}/genelists/Lee2016_RBFOX1.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/RBFOX1_targets.genes.list
#COSMIC gene lists
for class in oncogene all; do
  cat ${SFARI_ANNO}/genelists/COSMIC_census_${class}.genes.list | sort | uniq > \
  ${WRKDIR}/data/master_annotations/genelists/COSMIC_${class}.genes.list
done
cat ${SFARI_ANNO}/genelists/COSMIC_census_TSC.genes.list | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/COSMIC_tumor_suppressor.genes.list
#Collapse all tissue-specific gene sets to organ system-level sets
while read tissue; do
  echo ${tissue}
  #Iterate over geneset definitions
  while read geneset; do
    echo ${geneset}
    #Union
    while read file; do
      cat ${file}
    done < <( awk -v tissue=${tissue} -v geneset=${geneset} \
    '{ if ($1==tissue && $2==geneset) print $3 }' \
    ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list ) | \
    sort -k1,1 | uniq > \
    ${WRKDIR}/data/master_annotations/genelists/${tissue}_MASTER_UNION.${geneset}.genes.list
    #Intersection
    #Don't compute intersection for preferentially expressed genes
    # (By definition, preferential expression requires tissue-specificity per GTEx)
    if ! [ ${geneset} == "Preferentially_Expressed" ]; then
      if [ -e ${TMPDIR}/geneset_intersection.tmp ]; then
        rm ${TMPDIR}/geneset_intersection.tmp
      fi
      nlists=$( awk -v tissue=${tissue} -v geneset=${geneset} \
      '{ if ($1==tissue && $2==geneset) print $3 }' \
      ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list | wc -l )
      if [ ${nlists} -gt 0 ]; then
        while read file; do
          cat ${file}
        done < <( awk -v tissue=${tissue} -v geneset=${geneset} \
        '{ if ($1==tissue && $2==geneset) print $3 }' \
        ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list | head -n1 ) | \
        sort -k1,1 | uniq > \
        ${TMPDIR}/geneset_intersection.tmp
      fi
      if [ ${nlists} -gt 1 ]; then
        while read file; do
          fgrep -wf <( sed 's/\-/_/g' ${file} ) \
          <( sed 's/\-/_/g' ${TMPDIR}/geneset_intersection.tmp ) | sed 's/_/\-/g' > \
          ${TMPDIR}/geneset_intersection.tmp2
          mv ${TMPDIR}/geneset_intersection.tmp2 ${TMPDIR}/geneset_intersection.tmp
          wc -l ${TMPDIR}/geneset_intersection.tmp
        done < <( awk -v tissue=${tissue} -v geneset=${geneset} \
        '{ if ($1==tissue && $2==geneset) print $3 }' \
        ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list | sed '1d' )
      fi
      mv ${TMPDIR}/geneset_intersection.tmp \
      ${WRKDIR}/data/master_annotations/genelists/${tissue}_MASTER_INTERSECTION.${geneset}.genes.list
    fi
  done < <( awk -v tissue=${tissue} '{ if ($1==tissue) print $2 }' \
    ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list | \
    sort | uniq )
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_GeneSet_Linkers.list | \
  sort | uniq )
#Constrained AND highly expressed (at least one tissue)
cat ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | \
sort | uniq | sed 's/\-/_/g' | fgrep -wf \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | \
sort | uniq | sed 's/\_/\-/g' > \
${WRKDIR}/data/master_annotations/genelists/Constrained_AND_Highly_Expressed.genes.list
#Constrained OR highly expressed (at least one tissue)
cat ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | \
sort | uniq | sed 's/\-/_/g' | fgrep -wvf \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | \
sort | uniq | sed 's/\_/\-/g' > \
${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list
sort ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
uniq | sed 's/\-/_/g' | fgrep -wvf \
<( cat ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | \
   sed 's/\-/_/g' | sort | uniq ) | \
sort | uniq | sed 's/\_/\-/g' >> \
${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list
sort ${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list | \
uniq > ${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list2
mv ${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list2 \
${WRKDIR}/data/master_annotations/genelists/Constrained_OR_Highly_Expressed.genes.list
#Constrained NOT highly expressed (at least one tissue)
cat ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | \
sort | uniq | sed 's/\-/_/g' | fgrep -wvf \
<( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | \
sort | uniq | sed 's/\_/\-/g' > \
${WRKDIR}/data/master_annotations/genelists/Highly_Expressed_NOT_Constrained.genes.list
#Constrained NOT highly expressed (at least one tissue)
sort ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
uniq | sed 's/\-/_/g' | fgrep -wvf \
<( cat ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | \
   sed 's/\-/_/g' | sort | uniq ) | \
sort | uniq | sed 's/\_/\-/g' > \
${WRKDIR}/data/master_annotations/genelists/Constrained_NOT_Highly_Expressed.genes.list
#extTADA genes with q<0.05 from Nguyen et al., bioRxiv 2017 (Stahl/Sullivan labs)
#Note: Supplementary excel sheet downloaded from bioRxiv:
#  http://www.biorxiv.org/content/early/2017/08/25/135293
#Genes filtered on q<0.05 in excel for each phenotype and saved as gene lists on cluster
#DDG2P DD-associated genes from DECIPHER
#Note: must have manually downloaded csv, formatted it, and scp'ed it to
#  ERISOne at ${WRKDIR}/data/misc/DDG2P_genes.txt
awk -v FS="\t" '{ if (($3=="confirmed") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="dominant_negative" || $5=="loss_of_function")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Dominant_LOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="dominant_negative" || $5=="loss_of_function")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Dominant_LOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed") && \
                     ($4=="biallelic") && \
                     ($5=="dominant_negative" || $5=="loss_of_function")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Recessive_LOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
                     ($4=="biallelic") && \
                     ($5=="dominant_negative" || $5=="loss_of_function")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Recessive_LOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="activating" || $5=="increased_gene_dosage" || $5=="part_of_contiguous_gene_duplication")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Dominant_GOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="activating" || $5=="increased_gene_dosage" || $5=="part_of_contiguous_gene_duplication")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Dominant_GOF.genes.list
# awk -v FS="\t" '{ if (($3=="confirmed") && \
#                      ($4=="biallelic") && \
#                      ($5=="activating" || $5=="increased_gene_dosage" || $5=="part_of_contiguous_gene_duplication")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
# sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Recessive_GOF.genes.list
# awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
#                      ($4=="biallelic") && \
#                      ($5=="activating" || $5=="increased_gene_dosage" || $5=="part_of_contiguous_gene_duplication")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
# sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Recessive_GOF.genes.list
awk -v FS="\t" '{ if (($3=="confirmed") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="5_prime_or_3_prime_UTR_mutation" || $5=="" || $5=="all_missense/in_frame" || $5=="cis-regulatory_or_promotor_mutation" || $5=="uncertain")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Dominant_Unknown.genes.list
awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
                     ($4=="hemizygous" || $4=="monoallelic" || $4=="monoallelic,mosaic") && \
                     ($5=="5_prime_or_3_prime_UTR_mutation" || $5=="" || $5=="all_missense/in_frame" || $5=="cis-regulatory_or_promotor_mutation" || $5=="uncertain")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Dominant_Unknown.genes.list
awk -v FS="\t" '{ if (($3=="confirmed") && \
                     ($4=="biallelic") && \
                     ($5=="5_prime_or_3_prime_UTR_mutation" || $5=="" || $5=="all_missense/in_frame" || $5=="cis-regulatory_or_promotor_mutation" || $5=="uncertain")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Recessive_Unknown.genes.list
awk -v FS="\t" '{ if (($3=="confirmed" || $3=="probable" || $3=="possible") && \
                     ($4=="biallelic") && \
                     ($5=="5_prime_or_3_prime_UTR_mutation" || $5=="" || $5=="all_missense/in_frame" || $5=="cis-regulatory_or_promotor_mutation" || $5=="uncertain")) print $1 }' ${WRKDIR}/data/misc/DDG2P_genes.txt | \
sort | uniq > ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Recessive_Unknown.genes.list
#Coe 2017: recurrence of dnLGD and dnMIS from exomes on ~11k ASD/ID cases


#Get count of all genes and autosomal genes per gene list
while read list; do
  echo -e "${list}"
  #All genes
  cat ${list} | wc -l
  #Autosomal genes
  sed 's/\-/_/g' ${list} | fgrep -wf - \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) | \
  grep -e '^[0-9]' | cut -f4 | sort | uniq | wc -l
done < <( l ${WRKDIR}/data/master_annotations/genelists/*genes.list | \
  awk '{ print $9 }' | fgrep Coe2017 ) | paste - - -



##########################
##########################
#       NONCODING        #
##########################
##########################
#Conservative HARs from Kiana
#Note: must have uploaded from Slack into ${WRKDIR}/data/misc
sed '1d' ${WRKDIR}/data/misc/HARs.bed | cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > ${WRKDIR}/data/master_annotations/noncoding/HARs_conservative.elements.bed
#UCEs from Ruth (native hg19)
sed '1d' ${WRKDIR}/data/misc/hg19_Wren_UCEs_0based_full_set_5column_withheader.txt | \
cut -f1-3 | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/UCEs_McCole.elements.bed
#UCNEs from UCNEBase
cd ${WRKDIR}/data/misc/
wget ftp://ccg.vital-it.ch/UCNEbase/ucnes/hg19_UCNE_coord.bed
sed '1d' ${WRKDIR}/data/misc/hg19_UCNE_coord.bed | cut -f1-3 | sed 's/^chr//g' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/UCNEs.elements.bed
#Autosomal consensus syndromic CNV list
cat <( sed '1d' ${SFARI_ANNO}/misc/recurrent_CNV_loci_hg19.bed | sed 's/,//g' | cut -f1-3 ) \
<( sed '1d' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_CER.bed | cut -f1,4,5 ) \
${WRKDIR}/data/misc/PGC_CNV_all.bed | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
grep -e '^[0-9]' > ${WRKDIR}/data/master_annotations/noncoding/syndromic_CNVs.elements.bed
#FIREs
#Note: requires manually curated Schmitt_tissues.list
while read abbrev tissue; do
  fgrep -v "chrchr" ${SFARI_ANNO}/TADs/Schmitt2016/all_data_FIRE_calls/${abbrev}.FIRE.bed | \
  sed 's/^chr//g' > ${TMPDIR}/${abbrev}.FIRE.bed
  Rscript -e "options(scipen=1000); x <- read.table(\"${TMPDIR}/${abbrev}.FIRE.bed\",header=F); \
  write.table(x,\"${TMPDIR}/${abbrev}.FIRE.bed\",sep=\"\\t\",quote=F,col.names=F,row.names=F)"
  sort -Vk1,1 -k2,2n -k3,3n ${TMPDIR}/${abbrev}.FIRE.bed | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/FIREs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
#Merged FIREs
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/FIREs_${tissue}.elements.bed
done < <( fgrep PrimaryTissue ${WRKDIR}/data/misc/Schmitt_tissues.list ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d -1 -i - > \
${WRKDIR}/data/master_annotations/noncoding/FIREs_AllPrimaryTissues.elements.bed
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/FIREs_${tissue}.elements.bed
done < <( fgrep CellLine ${WRKDIR}/data/misc/Schmitt_tissues.list ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d -1 -i - > \
${WRKDIR}/data/master_annotations/noncoding/FIREs_AllCellLines.elements.bed
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/FIREs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d -1 -i - > \
${WRKDIR}/data/master_annotations/noncoding/FIREs_AllSources.elements.bed
#Conserved FIREs
awk -v OFS="\t" '{ print $1, $2, $3, 0 }' \
${WRKDIR}/data/master_annotations/noncoding/FIREs_AllSources.elements.bed > \
${TMPDIR}/FIREs_conserved.elements.bed
while read abbrev tissue; do
  bedtools intersect -c -a ${TMPDIR}/FIREs_conserved.elements.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/FIREs_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) print $1, $2, $3, $4+1; else print $1, $2, $3, $4 }' > \
  ${TMPDIR}/FIREs_conserved.elements.bed2
  mv ${TMPDIR}/FIREs_conserved.elements.bed2 ${TMPDIR}/FIREs_conserved.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
awk -v OFS="\t" '{ if ($4>=12) print $1, $2, $3 }' ${TMPDIR}/FIREs_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/FIREs_conserved.elements.bed
awk -v OFS="\t" '{ if ($4>=19) print $1, $2, $3 }' ${TMPDIR}/FIREs_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/FIREs_highlyConserved.elements.bed
#TBRs (Schmitt et al)
#Note: requires manually curated Schmitt_tissues.list
while read abbrev tissue; do
  sed 's/chr//g' ${SFARI_ANNO}/TADs/Schmitt2016/primary_cohort_TAD_boundaries/${abbrev}.IS.All_boundaries.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
#Merged TBRs
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed
done < <( fgrep PrimaryTissue ${WRKDIR}/data/misc/Schmitt_tissues.list ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 40000 -i - > \
${WRKDIR}/data/master_annotations/noncoding/TBRs_AllPrimaryTissues.elements.bed
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed
done < <( fgrep CellLine ${WRKDIR}/data/misc/Schmitt_tissues.list ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 40000 -i - > \
${WRKDIR}/data/master_annotations/noncoding/TBRs_AllCellLines.elements.bed
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 40000 -i - > \
${WRKDIR}/data/master_annotations/noncoding/TBRs_AllSources.elements.bed
#Conserved TBRs (Schmitt et al)
awk -v OFS="\t" '{ print $1, $2, $3, 0 }' \
${WRKDIR}/data/master_annotations/noncoding/TBRs_AllSources.elements.bed > \
${TMPDIR}/TBRs_conserved.elements.bed
while read abbrev tissue; do
  bedtools intersect -c -a ${TMPDIR}/TBRs_conserved.elements.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/TBRs_conserved.elements.bed2
  mv ${TMPDIR}/TBRs_conserved.elements.bed2 ${TMPDIR}/TBRs_conserved.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
awk -v OFS="\t" '{ if ($4>=12) print $1, $2, $3 }' ${TMPDIR}/TBRs_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/TBRs_conserved.elements.bed
awk -v OFS="\t" '{ if ($4>=19) print $1, $2, $3 }' ${TMPDIR}/TBRs_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/TBRs_highlyConserved.elements.bed
#TADs (Schmitt et al)
while read abbrev tissue; do
  paste ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed \
  <( sed '1d' ${WRKDIR}/data/master_annotations/noncoding/TBRs_${tissue}.elements.bed ) | \
  awk -v OFS="\t" '{ if ($1==$4) print $1, $3, $5 }' | awk '{ if ($3-$2<=5000000) print $0 }' > \
  ${WRKDIR}/data/master_annotations/noncoding/TADs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
#Compartments 
while read abbrev tissue; do
  echo ${tissue}
  if [ -e ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bw ]; then
    ${WRKDIR}/bin/utils/bigWigToBedGraph \
    ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bw \
    ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg;
  else
    lower=$( echo ${abbrev} | tr "A-Z" "a-z" )
    ${WRKDIR}/bin/utils/bigWigToBedGraph \
    ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${lower}.pc.bw \
    ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg
  fi
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
while read abbrev tissue; do
  cutoff=$( cut -f4 ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg | \
  sort -nrk1,1 | perl -e '$d=.1;@l=<>;print $l[int($d*$#l)]' )
  awk -v OFS="\t" -v cutoff=${cutoff} '{ if ($4>=cutoff) print $1, $2, $3 }' \
  ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | sed 's/^chr//g' > \
  ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsA_${tissue}.elements.bed
  cutoff=$( cut -f4 ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg | \
  sort -nrk1,1 | perl -e '$d=.9;@l=<>;print $l[int($d*$#l)]' )
  awk -v OFS="\t" -v cutoff=${cutoff} '{ if ($4<=cutoff) print $1, $2, $3 }' \
  ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${abbrev}.pc.bg | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | sed 's/^chr//g' > \
  ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsB_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
#Conserved A compartments
while read abbrev tissue; do
  cut -f1-3 ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsA_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - | awk -v OFS="\t" '{ print $0, "0" }' > ${TMPDIR}/strongCompartmentsA_conserved.elements.bed
while read abbrev tissue; do
  bedtools intersect -c -a ${TMPDIR}/strongCompartmentsA_conserved.elements.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsA_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/strongCompartmentsA_conserved.elements.bed2
  mv ${TMPDIR}/strongCompartmentsA_conserved.elements.bed2 \
  ${TMPDIR}/strongCompartmentsA_conserved.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
awk -v OFS="\t" '{ if ($4/21>0.5) print $1, $2, $3 }' ${TMPDIR}/strongCompartmentsA_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsA_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/21>=0.9) print $1, $2, $3 }' ${TMPDIR}/strongCompartmentsA_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsA_highlyConserved.elements.bed
#Conserved B compartments
while read abbrev tissue; do
  cut -f1-3 ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsB_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - | awk -v OFS="\t" '{ print $0, "0" }' > ${TMPDIR}/strongCompartmentsB_conserved.elements.bed
while read abbrev tissue; do
  bedtools intersect -c -a ${TMPDIR}/strongCompartmentsB_conserved.elements.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsB_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/strongCompartmentsB_conserved.elements.bed2
  mv ${TMPDIR}/strongCompartmentsB_conserved.elements.bed2 \
  ${TMPDIR}/strongCompartmentsB_conserved.elements.bed
done < ${WRKDIR}/data/misc/Schmitt_tissues.list
awk -v OFS="\t" '{ if ($4/21>0.5) print $1, $2, $3 }' ${TMPDIR}/strongCompartmentsB_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsB_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/21>=0.9) print $1, $2, $3 }' ${TMPDIR}/strongCompartmentsB_conserved.elements.bed > \
${WRKDIR}/data/master_annotations/noncoding/strongCompartmentsB_highlyConserved.elements.bed
#Enhancers (EnhancerAtlas)
mkdir ${WRKDIR}/data/misc/EnhancerAtlas
cd ${WRKDIR}/data/misc/EnhancerAtlas
wget http://enhanceratlas.org/data/enhseq/Astrocyte.fasta
wget http://enhanceratlas.org/data/enhseq/Bronchia_epithelial.fasta
wget http://enhanceratlas.org/data/enhseq/Esophagus.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_brain.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_heart.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_kidney.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_lung.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_muscle_leg.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_placenta.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_small_intestine.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_spinal_cord.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_stomach.fasta
wget http://enhanceratlas.org/data/enhseq/Fetal_thymus.fasta
wget http://enhanceratlas.org/data/enhseq/Fibroblast_foreskin.fasta
wget http://enhanceratlas.org/data/enhseq/Foreskin_keratinocyte.fasta
wget http://enhanceratlas.org/data/enhseq/Heart.fasta
wget http://enhanceratlas.org/data/enhseq/Left_ventricle.fasta
wget http://enhanceratlas.org/data/enhseq/Liver.fasta
wget http://enhanceratlas.org/data/enhseq/Lung.fasta
wget http://enhanceratlas.org/data/enhseq/Macrophage.fasta
wget http://enhanceratlas.org/data/enhseq/Myotube.fasta
wget http://enhanceratlas.org/data/enhseq/Osteoblast.fasta
wget http://enhanceratlas.org/data/enhseq/Ovary.fasta
wget http://enhanceratlas.org/data/enhseq/Pancreas.fasta
wget http://enhanceratlas.org/data/enhseq/Pancreatic_islet.fasta
wget http://enhanceratlas.org/data/enhseq/Skeletal_muscle.fasta
wget http://enhanceratlas.org/data/enhseq/Small_intestine.fasta
wget http://enhanceratlas.org/data/enhseq/Spleen.fasta
wget http://enhanceratlas.org/data/enhseq/Thymus.fasta
while read tissue; do
  echo "${tissue}"
  #CODE:
  # fgrep ">" ${WRKDIR}/data/misc/EnhancerAtlas/${tissue}.fasta | \
  # sed -e 's/>//g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/_/\t/g' -e 's/^chr//g' | \
  # awk -v OFS="\t" '{ if ($4=="") $4="."; print }' | \
  # sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - | \
  # sed -f ${WRKDIR}/data/master_annotations/gencode/ENSG_to_symbols.sed > \
  # ${WRKDIR}/data/master_annotations/noncoding/enhancers_${tissue}.elements.bed
  #PARALLELIZE:
  bsub -q short -sla miket_sc -u nobody -J ${tissue}_curateEnhancers \
  "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_enhancers.sh ${tissue}"
done < <( l ${WRKDIR}/data/misc/EnhancerAtlas/*fasta | awk '{ print $9 }' | \
  sed -e 's/\//\t/g' -e 's/\.fasta//g' | awk '{ print $NF }' )
#Conserved enhancers
while read tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/enhancers_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/EnhancerAtlas/*fasta | awk '{ print $9 }' | \
  sed -e 's/\//\t/g' -e 's/\.fasta//g' | awk '{ print $NF }' ) | cut -f1-3 | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" '{ print $1, $2, $3, 0 }' > \
${TMPDIR}/merged_enhancers.tmp
while read tissue; do
  bedtools intersect -c -a ${TMPDIR}/merged_enhancers.tmp \
  -b ${WRKDIR}/data/master_annotations/noncoding/enhancers_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/merged_enhancers.tmp2
  mv ${TMPDIR}/merged_enhancers.tmp2 ${TMPDIR}/merged_enhancers.tmp
done < <( l ${WRKDIR}/data/misc/EnhancerAtlas/*fasta | awk '{ print $9 }' | \
  sed -e 's/\//\t/g' -e 's/\.fasta//g' | awk '{ print $NF }' ) 
awk -v OFS="\t" '{ if ($4>15) print $1, $2, $3 }' ${TMPDIR}/merged_enhancers.tmp > \
${WRKDIR}/data/master_annotations/noncoding/enhancers_conserved.elements.bed
awk -v OFS="\t" '{ if ($4>26) print $1, $2, $3 }' ${TMPDIR}/merged_enhancers.tmp > \
${WRKDIR}/data/master_annotations/noncoding/enhancers_highlyConserved.elements.bed
#Super Enhancers (dbSUPER)
#Note: requires manually curated tissue file, ${WRKDIR}/data/misc/SuperEnhancer_tissues.list
mkdir ${WRKDIR}/data/misc/dbSUPER
cd ${WRKDIR}/data/misc/dbSUPER
wget http://bioinfo.au.tsinghua.edu.cn/dbsuper/data/bed/hg19/all_hg19_bed.zip
unzip ${WRKDIR}/data/misc/dbSUPER/all_hg19_bed.zip
while read tissue; do
  echo ${tissue}
  otissue=$( echo ${tissue} | sed 's/_/\ /g' )
  #Get super enhancers from file
  cut -f1-3 "${WRKDIR}/data/misc/dbSUPER/all_hg19_bed/${otissue}.bed" | \
  sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
  awk -v OFS="\t" '{ print $1, $2, $3, "SUPER_"NR }' > \
  ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed
  #Add all protein-coding genes within ±50kb
  awk -v OFS="\t" '{ print $1, $2-50000, $3+50000, $4 }' \
  ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($2<1) $2=1; print }' | bedtools intersect -wb -a - \
  -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.promoters.protein_coding.bed > \
  ${TMPDIR}/superEnh_wGenes.tmp
  while read chr start end enhID; do
    genes=$( fgrep -w ${enhID} ${TMPDIR}/superEnh_wGenes.tmp | cut -f9 | \
      sort | uniq | paste -s -d, | sed 's/\ //g' )
    if [ -z ${genes} ]; then
      genes="."
    fi
    echo -e "${chr}\t${start}\t${end}\t${genes}"
  done < ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed > \
  ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed2
  mv ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed2 \
  ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed
  #Rename file
  newname=$( echo "${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed" | \
  sed 's/\-//g' )
  if [ ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed != "${newname}" ]; then
    mv ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed \
    ${newname}
  fi
done < ${WRKDIR}/data/misc/SuperEnhancer_tissues.list
#Conserved super enhancers
while read tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed
done < <( sed 's/\-//g' ${WRKDIR}/data/misc/SuperEnhancer_tissues.list ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" '{ print $1, $2, $3, 0 }' > \
${TMPDIR}/merged_superEnhancers.tmp
while read tissue; do
  bedtools intersect -c -a ${TMPDIR}/merged_superEnhancers.tmp \
  -b ${WRKDIR}/data/master_annotations/noncoding/superEnhancers_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/merged_superEnhancers.tmp2
  mv ${TMPDIR}/merged_superEnhancers.tmp2 ${TMPDIR}/merged_superEnhancers.tmp
done < <( sed 's/\-//g' ${WRKDIR}/data/misc/SuperEnhancer_tissues.list ) 
awk -v OFS="\t" '{ if ($4>32) print $1, $2, $3 }' ${TMPDIR}/merged_superEnhancers.tmp > \
${WRKDIR}/data/master_annotations/noncoding/superEnhancers_conserved.elements.bed
awk -v OFS="\t" '{ if ($4>57) print $1, $2, $3 }' ${TMPDIR}/merged_superEnhancers.tmp > \
${WRKDIR}/data/master_annotations/noncoding/superEnhancers_highlyConserved.elements.bed
#ChromHMM on 75 epigenomes
#Note: requires manually curated tissue file, ${WRKDIR}/data/misc/Roadmap_Epi_tissues.list
#Note: also requires maunually curated ChromHMM mappings, 
mkdir ${WRKDIR}/data/misc/ChromHMM
cd ${WRKDIR}/data/misc/ChromHMM
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz
tar -xzvf ${WRKDIR}/data/misc/ChromHMM/all.mnemonics.bedFiles.tgz
while read EID tissue; do
  echo ${tissue}
  while read code state; do
    zcat ${WRKDIR}/data/misc/ChromHMM/${EID}_18_core_K27ac_mnemonics.bed.gz | \
    awk -v OFS="\t" -v code=${code} '{ if ($4==code) print $1, $2, $3 }' | \
    sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
    ${WRKDIR}/data/master_annotations/noncoding/${state}_${tissue}.elements.bed
  done < ${WRKDIR}/data/misc/ChromHMM_18way_states.list
done < ${WRKDIR}/data/misc/Roadmap_Epi_tissues.list
#Conserved ChromHMM for all states
while read code state; do
  echo ${state}
  while read EID tissue; do
    cat ${WRKDIR}/data/master_annotations/noncoding/${state}_${tissue}.elements.bed
  done < ${WRKDIR}/data/misc/Roadmap_Epi_tissues.list | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" \
  '{ print $1, $2, $3, 0 }' > ${TMPDIR}/${state}_merged.bed
  while read EID tissue; do
    bedtools intersect -c -a ${TMPDIR}/${state}_merged.bed \
    -b ${WRKDIR}/data/master_annotations/noncoding/${state}_${tissue}.elements.bed | \
    awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
    ${TMPDIR}/${state}_merged.bed2
    mv ${TMPDIR}/${state}_merged.bed2 ${TMPDIR}/${state}_merged.bed
    awk -v OFS="\t" '{ if ($4/75>0.5) print $1, $2, $3 }' ${TMPDIR}/${state}_merged.bed > \
    ${WRKDIR}/data/master_annotations/noncoding/${state}_conserved.elements.bed
    awk -v OFS="\t" '{ if ($4/75>0.9) print $1, $2, $3 }' ${TMPDIR}/${state}_merged.bed > \
    ${WRKDIR}/data/master_annotations/noncoding/${state}_highlyConserved.elements.bed
  done < ${WRKDIR}/data/misc/Roadmap_Epi_tissues.list
done < ${WRKDIR}/data/misc/ChromHMM_18way_states.list
#eQTLs by tissue (GTEx)
mkdir ${WRKDIR}/data/misc/GTEx_eQTL
cd ${WRKDIR}/data/misc/GTEx_eQTL
wget https://gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_eQTL.tar
tar -xvf ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL.tar
while read tissue; do
  #CODE:
  # ntissue=$( echo ${tissue} | sed 's/\-/_/g' )
  # zcat ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL/${tissue}_Analysis.v6p.signif_snpgene_pairs.txt.gz | \
  # cut -f1-2 | sed 's/\./\t/g' | cut -f1-2 | sed '1d' | sed 's/_/\t/g' | awk -v OFS="\t" \
  # '{ print $1, $2, $2+length($4), $6 }' | sort -Vk1,1 -k2,2n -k3,3n | \
  # bedtools merge -c 4 -o distinct -i - | \
  # sed -f ${WRKDIR}/data/master_annotations/gencode/ENSG_to_symbols.sed > \
  # ${WRKDIR}/data/master_annotations/noncoding/eQTLs_${ntissue}.elements.bed
  #PARALLELIZE:
  bsub -q short -sla miket_sc -u nobody -J ${tissue}_eQTLs \
  "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_eQTLs.sh ${tissue}"
done < <( l ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL/*_Analysis.v6p.signif_snpgene_pairs.txt.gz | \
  sed 's/\//\t/g' | awk '{ print $NF }' | \
  sed 's/_Analysis\.v6p\.signif_snpgene_pairs\.txt\.gz//g' | sort -Vk1,1 )
#Conserved eQTLs
while read tissue; do
  ntissue=$( echo ${tissue} | sed 's/\-/_/g' )
  cut -f1-3 ${WRKDIR}/data/master_annotations/noncoding/eQTLs_${ntissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL/*_Analysis.v6p.signif_snpgene_pairs.txt.gz | \
  sed 's/\//\t/g' | awk '{ print $NF }' | \
  sed 's/_Analysis\.v6p\.signif_snpgene_pairs\.txt\.gz//g' | sort -Vk1,1 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" \
'{ print $1, $2, $3, "0" }' > ${TMPDIR}/eQTLs_merged.bed
while read tissue; do
  ntissue=$( echo ${tissue} | sed 's/\-/_/g' )
  bedtools intersect -c -a ${TMPDIR}/eQTLs_merged.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/eQTLs_${ntissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/eQTLs_merged.bed2
  mv ${TMPDIR}/eQTLs_merged.bed2 ${TMPDIR}/eQTLs_merged.bed
done < <( l ${WRKDIR}/data/misc/GTEx_eQTL/GTEx_Analysis_v6p_eQTL/*_Analysis.v6p.signif_snpgene_pairs.txt.gz | \
  sed 's/\//\t/g' | awk '{ print $NF }' | \
  sed 's/_Analysis\.v6p\.signif_snpgene_pairs\.txt\.gz//g' | sort -Vk1,1 )
awk -v OFS="\t" '{ if ($4/44>0.5) print $1, $2, $3 }' ${TMPDIR}/eQTLs_merged.bed | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/eQTLs_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/44>0.9) print $1, $2, $3 }' ${TMPDIR}/eQTLs_merged.bed | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/eQTLs_highlyConserved.elements.bed
#CpG Islands
cd ${WRKDIR}/data/misc
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExtUnmasked.txt.gz
zcat ${WRKDIR}/data/misc/cpgIslandExtUnmasked.txt.gz | cut -f2-4 | sed 's/^chr//g' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/CpG_islands.elements.bed
#PhyloP conservation peaks
mkdir ${WRKDIR}/data/misc/PhyloP
cd ${WRKDIR}/data/misc/PhyloP
for chr in $( seq 1 22 ); do
  bsub -q filemove -sla miket_sc -J PhyloP_DL -u nobody \
  "wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr${chr}.phyloP46way.placental.wigFix.gz `pwd`"
done
for chr in $( seq 1 22 ); do
  ${WRKDIR}/bin/utils/wigToBigWig \
  ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.wigFix.gz \
  <( sed 's/^/chr/g' /data/talkowski/rlc47/src/GRCh37.genome ) \
  ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bigWig
done 
for chr in $( seq 1 22 ); do
  echo ${chr}
  ${WRKDIR}//bin/utils/bigWigToBedGraph -chrom=chr${chr} \
  ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bigWig \
  ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bg
done 
for chr in $( seq 1 22 ); do
  #CODE:
  # awk -v OFS="\t" '{ if ($4>=1) print $1, $2, $3 }' \
  # ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bg | \
  # sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 100 -i - | \
  # awk -v OFS="\t" '{ if ($3-$2>=200) print $1, $2, $3 }' | sed 's/^chr//g'
  #PARALLELIZE:
  bsub -q short -sla miket_sc -u nobody -J PhyloP_chr${chr} \
  "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_PhyloP.sh ${chr}"
done 
for chr in $( seq 1 22 ); do
  cat ${WRKDIR}/data/misc/PhyloP/evolutionarilyConserved_PhyloP.chr${chr}.elements.bed
done | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/master_annotations/noncoding/evolutionarilyConserved_PhyloP.elements.bed
#GERP conservation peaks
mkdir ${WRKDIR}/data/misc/GERP
cd ${WRKDIR}/data/misc/GERP
wget http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_elements.tar.gz
tar -xzvf ${WRKDIR}/data/misc/hg19.GERP_elements.tar.gz
for chr in $( seq 1 22 ); do
  cat ${WRKDIR}/data/misc/GERP/hg19_chr${chr}_elems.txt | \
  awk -v OFS="\t" -v chr=${chr} '{ print chr, $1, $2 }'
done | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/evolutionarilyConserved_GERP.elements.bed
#Consensus conservation peaks
bedtools intersect -r -f 0.1 -wa -wb \
-a ${WRKDIR}/data/master_annotations/noncoding/evolutionarilyConserved_PhyloP.elements.bed \
-b ${WRKDIR}/data/master_annotations/noncoding/evolutionarilyConserved_GERP.elements.bed | \
awk -v OFS="\t" '{ print $1, $2, $3"\n"$4, $5, $6 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/evolutionarilyConserved_consensus.elements.bed
#Replication timing
cd ${WRKDIR}/data/misc
wget http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip
unzip ${WRKDIR}/data/misc/Koren-et-al-Table-S2.zip
mv ${WRKDIR}/data/misc/Koren\ et\ al\ Table\ S2.txt \
${WRKDIR}/data/misc/Koren_et_al_Table_S2.txt
paste ${WRKDIR}/data/misc/Koren_et_al_Table_S2.txt \
<( sed '1d' ${WRKDIR}/data/misc/Koren_et_al_Table_S2.txt ) | awk -v OFS="\t" \
'{ if ($1!=$4) $2=0; print $4, $2, $5, $6 }' | awk -v OFS="\t" \
'{ if ($4>=1 && $4!="NaN" && $3>=$2) printf "%i\t%.0f\t%.0f\n", $1, $2, $3 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -d 10000 -i - > \
${WRKDIR}/data/master_annotations/noncoding/earlyReplicating_Koren.elements.bed
paste ${WRKDIR}/data/misc/Koren_et_al_Table_S2.txt \
<( sed '1d' ${WRKDIR}/data/misc/Koren_et_al_Table_S2.txt ) | awk -v OFS="\t" \
'{ if ($1!=$4) $2=0; print $4, $2, $5, $6 }' | awk -v OFS="\t" \
'{ if ($4<=-1 && $4!="NaN" && $3>=$2) printf "%i\t%.0f\t%.0f\n", $1, $2, $3 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -d 10000 -i - > \
${WRKDIR}/data/master_annotations/noncoding/lateReplicating_Koren.elements.bed
#Differentially methylated regions (DMRs)
mkdir ${WRKDIR}/data/misc/DMRs
cd ${WRKDIR}/data/misc/DMRs
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_SamplesOverview.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_SupplementaryTable.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_SupplementaryTable_readme.txt.gz
${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_DMRs.R \
${WRKDIR}/data/misc/DMRs/GSE46644_SupplementaryTable.txt.gz \
${WRKDIR}/data/misc/DMRs/
while read tissue; do
  sort -Vk1,1 -k2,2n -k3,3n ${WRKDIR}/data/misc/DMRs/${tissue}_mean.DMRs.bed | \
  bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/DMRs_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/DMRs/*bed | awk '{ print $9 }' | cut -f1 -d\. | \
  sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/_mean//g' )
#Partially methylated domains (PMDs)
#Note: had to manually curate PMDs per tissue from Supplementary Table 8 of
#Schultz et al., Nature (2015)
#Also needs manually curated abbreviation-to-tissue linking file
while read abbrev tissue; do
  echo ${tissue}
  cat ${WRKDIR}/data/misc/PMDs_curated/${abbrev}*bed | \
  sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/master_annotations/noncoding/PMDs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schultz_PMD_tissues.list
#Conserved PMDs
while read abbrev tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/PMDs_${tissue}.elements.bed
done < ${WRKDIR}/data/misc/Schultz_PMD_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - | awk -v OFS="\t" '{ print $1, $2, $3, "0" }' > \
${TMPDIR}/PMDs_conserved.tmp
while read abbrev tissue; do
  bedtools intersect -c -a ${TMPDIR}/PMDs_conserved.tmp \
  -b ${WRKDIR}/data/master_annotations/noncoding/PMDs_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/PMDs_conserved.tmp2
  mv ${TMPDIR}/PMDs_conserved.tmp2 ${TMPDIR}/PMDs_conserved.tmp
done < ${WRKDIR}/data/misc/Schultz_PMD_tissues.list
awk -v OFS="\t" '{ if ($4/18>0.5) print $1, $2, $3 }' ${TMPDIR}/PMDs_conserved.tmp > \
${WRKDIR}/data/master_annotations/noncoding/PMDs_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/18>0.9) print $1, $2, $3 }' ${TMPDIR}/PMDs_conserved.tmp > \
${WRKDIR}/data/master_annotations/noncoding/PMDs_highlyConserved.elements.bed
#GWAS catalogue (split by disease phenotype)
#Note: must have downloaded GWAS catalogue and moved it to cluster already
#Must have also cleaned GWAS catalogue and moved phenotype matrix on cluster
sed 's/\ /_/g' ${WRKDIR}/data/misc/gwas_catalog_v1.0.1-associations_e88_r2017-04-03.tsv | \
awk -v FS="\t" -v OFS="\t" '{ if ($34=="N" && $29>=8) print $0 }' | awk -v FS="\t" -v OFS="\t" \
'$1 ~ /2013|2014|2015|2016|2017/ { print $12, $13, $35 }' | \
awk -v OFS="\t" -v FS="\t" '$1 !~ /x/ { print $1, $2, $3 }' | \
awk -v OFS="\t" -v FS="\t" '{ if ($1!="" && $2!="" && $3!="") print $0 }' > \
${TMPDIR}/GWAS_catalogue.formatted.tmp
awk -v OFS="\t" -v FS="\t" '$1 !~ /\;/ { print $1, $2, $2+1, $3 }' \
${TMPDIR}/GWAS_catalogue.formatted.tmp > \
${WRKDIR}/data/misc/gwas_catalog_cleaned.sites_and_phenos.bed
while read chrs starts skip pheno; do
  chr=$( echo ${chrs} | sed 's/\;/\n/g' | sort | uniq | head -n1 )
  echo "${starts}" | sed 's/\;/\n/g' | awk -v OFS="\t" -v chr=${chr} -v pheno=${pheno} \
  '{ print chr, $1, $1+2, pheno }'
done < <( awk -v OFS="\t" -v FS="\t" '$1 ~ /\;/ { print $1, $2, $2+1, $3 }' \
  ${TMPDIR}/GWAS_catalogue.formatted.tmp ) >> \
${WRKDIR}/data/misc/gwas_catalog_cleaned.sites_and_phenos.bed
while read pheno; do
  echo ${pheno}
  #Get column according to phenotype
  col=$( head -n1 ${WRKDIR}/data/misc/GWAS_catalogue_phenoInclusionMatrix.txt | \
    sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w ${pheno} | cut -f1 )
  #Get terms associated with that column and pull associated sites, add ±10kb, and merge
  awk -v col=${col} -v OFS="\t" '{ if ($(col)==1) print $1 }' \
  ${WRKDIR}/data/misc/GWAS_catalogue_phenoInclusionMatrix.txt | \
  fgrep -f - ${WRKDIR}/data/misc/gwas_catalog_cleaned.sites_and_phenos.bed | \
  awk -v OFS="\t" '{ print "chr"$1, $2-5000, $3+5000 }' | awk -v OFS="\t" \
  '{ if ($2<1) $2=1; print }' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${TMPDIR}/${pheno}.loci.hg38.bed
  #Lift over from hg38 to GRCh37
  liftOver -minMatch=0.5 ${TMPDIR}/${pheno}.loci.hg38.bed \
  /data/talkowski/rlc47/src/hg38ToHg19.over.chain.gz ${TMPDIR}/${pheno}.loci.hg19.bed \
  ${TMPDIR}/${pheno}.loci.hg38_liftOverFail.bed
  #Print stats on loci that failed remapping
  cat ${TMPDIR}/${pheno}.loci.hg19.bed | wc -l
  cat ${TMPDIR}/${pheno}.loci.hg38.bed | wc -l
  #Clean and write to file
  sed 's/^chr//g' ${TMPDIR}/${pheno}.loci.hg19.bed | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/master_annotations/noncoding/GWAS_loci_${pheno}.elements.bed
done < <( head -n1 ${WRKDIR}/data/misc/GWAS_catalogue_phenoInclusionMatrix.txt | \
  cut -f2- | sed 's/\t/\n/g' ) | fgrep -v Reading | fgrep -v Mapping | paste - - -
#Segmental duplications
cd ${WRKDIR}/data/misc/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/genomicSuperDups.txt.gz
zcat ${WRKDIR}/data/misc/genomicSuperDups.txt.gz | awk -v OFS="\t" \
'{ print $2, $3, $4 }' | sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/SegDups.elements.bed
#Repeat masker
cd ${WRKDIR}/data/misc/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz
while read class; do
  echo ${class}
  zcat ${WRKDIR}/data/misc/rmsk.txt.gz | awk -v OFS="\t" -v class=${class} \
  '{ if ($12==class) print $6, $7, $8 }' | sed 's/^chr//g' | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/RepeatMasker_${class}.elements.bed
done < <( zcat ${WRKDIR}/data/misc/rmsk.txt.gz | cut -f12 | sort | uniq | \
  grep -ve 'Other\|RNA\|Unknown\|\?' )
#Recombination frequency
cd ${WRKDIR}/data/misc/
wget https://www.decode.com/additional/sex-averaged.rmap
cutoff=$( sed '1d' ${WRKDIR}/data/misc/sex-averaged.rmap | awk '{ if ($3>0) print $4 }' | \
sort -nrk1,1 | perl -e '$d=.1;@l=<>;print $l[int($d*$#l)]' | sed 's/\ //g' )
sed '1d' ${WRKDIR}/data/misc/sex-averaged.rmap | awk -v OFS="\t" -v cutoff=${cutoff} \
'{ if ($3>0 && $4>=cutoff) print $1, $2-5000, $2+5000 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > ${TMPDIR}/recomb_hotspots.hg18.bed
liftOver -minMatch=0.5 ${TMPDIR}/recomb_hotspots.hg18.bed \
/data/talkowski/rlc47/src/hg18ToHg19.over.chain ${TMPDIR}/recomb_hotspots.hg19.bed \
${TMPDIR}/recomb_hotspots.hg18_liftOverFail.bed
sed 's/^chr//g' ${TMPDIR}/recomb_hotspots.hg19.bed | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > \
${WRKDIR}/data/master_annotations/noncoding/Recombination_Hotspots.elements.bed
#DHS - primary tissues & primary cells
mkdir ${WRKDIR}/data/misc/DHS
cd ${WRKDIR}/data/misc/DHS
wget -O ./FetalArmMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF852IEJ/@@download/ENCFF852IEJ.bed.gz
wget -O ./FetalHeart_DHS.bed.gz https://www.encodeproject.org/files/ENCFF559XLA/@@download/ENCFF559XLA.bed.gz
wget -O ./Heart_DHS.bed.gz https://www.encodeproject.org/files/ENCFF822VOH/@@download/ENCFF822VOH.bed.gz
wget -O ./FetalStomach_DHS.bed.gz https://www.encodeproject.org/files/ENCFF191JBX/@@download/ENCFF191JBX.bed.gz
wget -O ./Stomach_DHS.bed.gz https://www.encodeproject.org/files/ENCFF701AFB/@@download/ENCFF701AFB.bed.gz
wget -O ./FetalLargeIntestine_DHS.bed.gz https://www.encodeproject.org/files/ENCFF539JKS/@@download/ENCFF539JKS.bed.gz
wget -O ./FetalBackMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF317QQX/@@download/ENCFF317QQX.bed.gz
wget -O ./FetalLegMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF332HJI/@@download/ENCFF332HJI.bed.gz
wget -O ./FetalLeftLung_DHS.bed.gz https://www.encodeproject.org/files/ENCFF473MAD/@@download/ENCFF473MAD.bed.gz
wget -O ./FetalSmallIntestine_DHS.bed.gz https://www.encodeproject.org/files/ENCFF866PBJ/@@download/ENCFF866PBJ.bed.gz
wget -O ./FetalBrain_DHS.bed.gz https://www.encodeproject.org/files/ENCFF251UGS/@@download/ENCFF251UGS.bed.gz
wget -O ./FetalKidney_DHS.bed.gz https://www.encodeproject.org/files/ENCFF475AOQ/@@download/ENCFF475AOQ.bed.gz
wget -O ./FetalLung_DHS.bed.gz https://www.encodeproject.org/files/ENCFF605MPK/@@download/ENCFF605MPK.bed.gz
wget -O ./FetalRenalCortexInterstitium_DHS.bed.gz https://www.encodeproject.org/files/ENCFF617JOF/@@download/ENCFF617JOF.bed.gz
wget -O ./RenalCortexInterstitium_DHS.bed.gz https://www.encodeproject.org/files/ENCFF914YLY/@@download/ENCFF914YLY.bed.gz
wget -O ./FetalRightLung_DHS.bed.gz https://www.encodeproject.org/files/ENCFF828TWL/@@download/ENCFF828TWL.bed.gz
wget -O ./AdrenalGland_DHS.bed.gz https://www.encodeproject.org/files/ENCFF488DIV/@@download/ENCFF488DIV.bed.gz
wget -O ./FetalAdrenalGland_DHS.bed.gz https://www.encodeproject.org/files/ENCFF189FNX/@@download/ENCFF189FNX.bed.gz
wget -O ./RenalPelvis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF958VNA/@@download/ENCFF958VNA.bed.gz
wget -O ./FetalRenalPelvis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF306IHB/@@download/ENCFF306IHB.bed.gz
wget -O ./FetalRightKidney_DHS.bed.gz https://www.encodeproject.org/files/ENCFF593JZH/@@download/ENCFF593JZH.bed.gz
wget -O ./FetalThymus_DHS.bed.gz https://www.encodeproject.org/files/ENCFF196JFB/@@download/ENCFF196JFB.bed.gz
wget -O ./FetalLeftKidney_DHS.bed.gz https://www.encodeproject.org/files/ENCFF246DEW/@@download/ENCFF246DEW.bed.gz
wget -O ./Placenta_DHS.bed.gz https://www.encodeproject.org/files/ENCFF012LUH/@@download/ENCFF012LUH.bed.gz
wget -O ./FetalSpinalCord_DHS.bed.gz https://www.encodeproject.org/files/ENCFF795ISF/@@download/ENCFF795ISF.bed.gz
wget -O ./FetalLeftRenalCortexInterstitium_DHS.bed.gz https://www.encodeproject.org/files/ENCFF622RNC/@@download/ENCFF622RNC.bed.gz
wget -O ./FetalLeftRenalPelvis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF802SNU/@@download/ENCFF802SNU.bed.gz
wget -O ./FetalRightRenalPelvis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF242ABW/@@download/ENCFF242ABW.bed.gz
wget -O ./FetalTrunkMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF484TJS/@@download/ENCFF484TJS.bed.gz
wget -O ./Ovary_DHS.bed.gz https://www.encodeproject.org/files/ENCFF119HQB/@@download/ENCFF119HQB.bed.gz
wget -O ./FetalRetina_DHS.bed.gz https://www.encodeproject.org/files/ENCFF576YYS/@@download/ENCFF576YYS.bed.gz
wget -O ./FetalRightRenalCortexInterstitium_DHS.bed.gz https://www.encodeproject.org/files/ENCFF408KGV/@@download/ENCFF408KGV.bed.gz
wget -O ./Testis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF785AEY/@@download/ENCFF785AEY.bed.gz
wget -O ./FetalTestis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF502RSX/@@download/ENCFF502RSX.bed.gz
wget -O ./Thyroid_DHS.bed.gz https://www.encodeproject.org/files/ENCFF563DYP/@@download/ENCFF563DYP.bed.gz
wget -O ./TransverseColon_DHS.bed.gz https://www.encodeproject.org/files/ENCFF697ZFE/@@download/ENCFF697ZFE.bed.gz
wget -O ./PancreasBody_DHS.bed.gz https://www.encodeproject.org/files/ENCFF454JYJ/@@download/ENCFF454JYJ.bed.gz
wget -O ./FetalFacialProminence_DHS.bed.gz https://www.encodeproject.org/files/ENCFF243KHX/@@download/ENCFF243KHX.bed.gz
wget -O ./FetalEye_DHS.bed.gz https://www.encodeproject.org/files/ENCFF504ZVA/@@download/ENCFF504ZVA.bed.gz
wget -O ./FetalLimb_DHS.bed.gz https://www.encodeproject.org/files/ENCFF423FPL/@@download/ENCFF423FPL.bed.gz
wget -O ./Pancreas_DHS.bed.gz https://www.encodeproject.org/files/ENCFF890KUA/@@download/ENCFF890KUA.bed.gz
wget -O ./SigmoidColon_DHS.bed.gz https://www.encodeproject.org/files/ENCFF693FUF/@@download/ENCFF693FUF.bed.gz
wget -O ./FetalTongue_DHS.bed.gz https://www.encodeproject.org/files/ENCFF831SBJ/@@download/ENCFF831SBJ.bed.gz
wget -O ./Uterus_DHS.bed.gz https://www.encodeproject.org/files/ENCFF645GPP/@@download/ENCFF645GPP.bed.gz
wget -O ./AmmonsHorn_DHS.bed.gz https://www.encodeproject.org/files/ENCFF236UMJ/@@download/ENCFF236UMJ.bed.gz
wget -O ./CerebellarCortex_DHS.bed.gz https://www.encodeproject.org/files/ENCFF482VEA/@@download/ENCFF482VEA.bed.gz
wget -O ./EsophagusSquamousEpithelium_DHS.bed.gz https://www.encodeproject.org/files/ENCFF570SFB/@@download/ENCFF570SFB.bed.gz
wget -O ./FetalForelimb_DHS.bed https://www.encodeproject.org/files/ENCFF268PYO/@@download/ENCFF268PYO.bed.gz
wget -O ./GastrocnemiusMedialis_DHS.bed.gz https://www.encodeproject.org/files/ENCFF067SRV/@@download/ENCFF067SRV.bed.gz
wget -O ./GlobusPallidus_DHS.bed.gz https://www.encodeproject.org/files/ENCFF461EHY/@@download/ENCFF461EHY.bed.gz
wget -O ./FetalLeftVentricle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF315TJF/@@download/ENCFF315TJF.bed.gz
wget -O ./FetalHindlimbMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF745LBF/@@download/ENCFF745LBF.bed.gz
wget -O ./InferiorParietalCortex_DHS.bed.gz https://www.encodeproject.org/files/ENCFF913FRG/@@download/ENCFF913FRG.bed.gz
wget -O ./IsletOfLangerhans_DHS.bed.gz https://www.encodeproject.org/files/ENCFF001UXL/@@download/ENCFF001UXL.bed.gz
wget -O ./MedullaOblongata_DHS.bed.gz https://www.encodeproject.org/files/ENCFF177LYE/@@download/ENCFF177LYE.bed.gz
wget -O ./Midbrain_DHS.bed.gz https://www.encodeproject.org/files/ENCFF490EZU/@@download/ENCFF490EZU.bed.gz
wget -O ./MiddleFrontalGyrus_DHS.bed.gz https://www.encodeproject.org/files/ENCFF440HVH/@@download/ENCFF440HVH.bed.gz
wget -O ./OccipitalLobe_DHS.bed.gz https://www.encodeproject.org/files/ENCFF261WWC/@@download/ENCFF261WWC.bed.gz
wget -O ./Pons_DHS.bed.gz https://www.encodeproject.org/files/ENCFF615GGO/@@download/ENCFF615GGO.bed.gz
wget -O ./Prostate_DHS.bed.gz https://www.encodeproject.org/files/ENCFF267AWH/@@download/ENCFF267AWH.bed.gz
wget -O ./LiverRightLobe_DHS.bed.gz https://www.encodeproject.org/files/ENCFF865CGA/@@download/ENCFF865CGA.bed.gz
wget -O ./Skin_DHS.bed.gz https://www.encodeproject.org/files/ENCFF976WFS/@@download/ENCFF976WFS.bed.gz
wget -O ./Spleen_DHS.bed.gz https://www.encodeproject.org/files/ENCFF382VKO/@@download/ENCFF382VKO.bed.gz
wget -O ./SuperiorTemporalGyrus_DHS.bed.gz https://www.encodeproject.org/files/ENCFF835REN/@@download/ENCFF835REN.bed.gz
wget -O ./FetalThoracicMuscle_DHS.bed.gz https://www.encodeproject.org/files/ENCFF917UMO/@@download/ENCFF917UMO.bed.gz
wget -O ./UmbilicalCord_DHS.bed.gz https://www.encodeproject.org/files/ENCFF442KIW/@@download/ENCFF442KIW.bed.gz
wget -O ./FetalUrinaryBladder_DHS.bed.gz https://www.encodeproject.org/files/ENCFF230RHI/@@download/ENCFF230RHI.bed.gz
wget -O ./Bcell_DHS.bed.gz https://www.encodeproject.org/files/ENCFF001VYW/@@download/ENCFF001VYW.bed.gz
wget -O ./Tcell_DHS.bed.gz https://www.encodeproject.org/files/ENCFF026SFK/@@download/ENCFF026SFK.bed.gz
wget -O ./Monocyte_DHS.bed.gz https://www.encodeproject.org/files/ENCFF564MST/@@download/ENCFF564MST.bed.gz
wget -O ./NKcell_DHS.bed.gz https://www.encodeproject.org/files/ENCFF224LJW/@@download/ENCFF224LJW.bed.gz
zcat ${WRKDIR}/data/misc/DHS/SuperiorTemporalGyrus_DHS.bed.gz | awk -v OFS="\t" \
'{ print $1, $2, $3, $4, $5, $6, $9, $7, $8, $9 }' > \
${WRKDIR}/data/misc/DHS/SuperiorTemporalGyrus_DHS.bed
gzip -f ${WRKDIR}/data/misc/DHS/SuperiorTemporalGyrus_DHS.bed
${WRKDIR}/data/misc/DHS/SuperiorTemporalGyrus_DHS.bed
while read tissue; do
  echo ${tissue}
  #All DHS
  zcat ${WRKDIR}/data/misc/DHS/${tissue}_DHS.bed.gz | sed 's/^chr//g' | \
  cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/DHS_${tissue}.elements.bed
  #Strong DHS (top 10%)
  cutoff=$( zcat ${WRKDIR}/data/misc/DHS/${tissue}_DHS.bed.gz | cut -f7 | \
  sort -nrk1,1 | perl -e '$d=.1;@l=<>;print $l[int($d*$#l)]' )
  zcat ${WRKDIR}/data/misc/DHS/${tissue}_DHS.bed.gz | awk -v OFS="\t" -v cutoff=${cutoff} \
  '{ if ($7>=cutoff) print $1, $2, $3 }' | sed 's/^chr//g' | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/StrongDHS_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/DHS/*_DHS.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_DHS/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' )
#Conserved DHS
while read tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/DHS_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/DHS/*_DHS.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_DHS/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" \
  '{ print $1, $2, $3, "0" }' > ${TMPDIR}/DHS_merged.bed
while read tissue; do
  bedtools intersect -c -a ${TMPDIR}/DHS_merged.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/DHS_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/DHS_merged.bed2
  mv ${TMPDIR}/DHS_merged.bed2 ${TMPDIR}/DHS_merged.bed
done < <( l ${WRKDIR}/data/misc/DHS/*_DHS.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_DHS/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' )
awk -v OFS="\t" '{ if ($4/72>0.5) print $1, $2, $3 }' ${TMPDIR}/DHS_merged.bed > \
${WRKDIR}/data/master_annotations/noncoding/DHS_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/72>0.9) print $1, $2, $3 }' ${TMPDIR}/DHS_merged.bed > \
${WRKDIR}/data/master_annotations/noncoding/DHS_highlyConserved.elements.bed
#H3K27ac peaks
mkdir ${WRKDIR}/data/misc/H3K27ac/
cd ${WRKDIR}/data/misc/H3K27ac/
wget -O ./FetalAdrenalGland_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF433MRW/@@download/ENCFF433MRW.bed.gz
wget -O ./AdrenalGland_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF812UKC/@@download/ENCFF812UKC.bed.gz
wget -O ./Pancreas_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF821IQZ/@@download/ENCFF821IQZ.bed.gz
wget -O ./GastrocnemiusMedialis_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF090ODW/@@download/ENCFF090ODW.bed.gz
wget -O ./FetalSmallIntestine_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF264TVA/@@download/ENCFF264TVA.bed.gz
wget -O ./SmallIntestine_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF452ODQ/@@download/ENCFF452ODQ.bed.gz
wget -O ./Spleen_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF434AUI/@@download/ENCFF434AUI.bed.gz
wget -O ./FetalStomach_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF720RQS/@@download/ENCFF720RQS.bed.gz
wget -O ./Stomach_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF868AQH/@@download/ENCFF868AQH.bed.gz
wget -O ./Thyroid_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF093DSS/@@download/ENCFF093DSS.bed.gz
wget -O ./LeftVentricle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF546TJN/@@download/ENCFF546TJN.bed.gz
wget -O ./RectalMucosa_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF030XSU/@@download/ENCFF030XSU.bed.gz
wget -O ./Psoas_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF761BZK/@@download/ENCFF761BZK.bed.gz
wget -O ./ThoracicAorta_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF737DEL/@@download/ENCFF737DEL.bed.gz
wget -O ./FetalThymus_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF419EUC/@@download/ENCFF419EUC.bed.gz
wget -O ./PeyersPatch_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF357MJZ/@@download/ENCFF357MJZ.bed.gz
wget -O ./Adipose_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF222CHB/@@download/ENCFF222CHB.bed.gz
wget -O ./FetalAmnion_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF206LAP/@@download/ENCFF206LAP.bed.gz
wget -O ./BreastEpithelium_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF154XFN/@@download/ENCFF154XFN.bed.gz
wget -O ./ColonicMucosa_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF723NUR/@@download/ENCFF723NUR.bed.gz
wget -O ./CoronaryArtery_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF130AWQ/@@download/ENCFF130AWQ.bed.gz
wget -O ./EndocrinePancreas_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF885UNK/@@download/ENCFF885UNK.bed.gz
wget -O ./Esophagus_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF694EXU/@@download/ENCFF694EXU.bed.gz
wget -O ./EsophagusMuscularisMucosa_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF903OJB/@@download/ENCFF903OJB.bed.gz
wget -O ./EsophagusSquamousEpithelium_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF187KWD/@@download/ENCFF187KWD.bed.gz
wget -O ./GastroesophagealSphincter_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF177WNT/@@download/ENCFF177WNT.bed.gz
wget -O ./RightVentricle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF995GJW/@@download/ENCFF995GJW.bed.gz
wget -O ./FetalLargeIntestine_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF906VKA/@@download/ENCFF906VKA.bed.gz
wget -O ./Lung_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF768OEM/@@download/ENCFF768OEM.bed.gz
wget -O ./ColonMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF740FVQ/@@download/ENCFF740FVQ.bed.gz
wget -O ./DuodenumMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF063BDL/@@download/ENCFF063BDL.bed.gz
wget -O ./FetalLegMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF563PBB/@@download/ENCFF563PBB.bed.gz
wget -O ./FetalTrunkMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF711AVI/@@download/ENCFF711AVI.bed.gz
wget -O ./Ovary_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF546HYT/@@download/ENCFF546HYT.bed.gz
wget -O ./FetalPlacenta_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF107GTA/@@download/ENCFF107GTA.bed.gz
wget -O ./RectalSmoothMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF123SVJ/@@download/ENCFF123SVJ.bed.gz
wget -O ./RightAtriumAuricularRegion_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF159NYF/@@download/ENCFF159NYF.bed.gz
wget -O ./RightAtrium_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF101VMB/@@download/ENCFF101VMB.bed.gz
wget -O ./LiverRightLobe_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF285MWK/@@download/ENCFF285MWK.bed.gz
wget -O ./SigmoidColon_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF027RPB/@@download/ENCFF027RPB.bed.gz
wget -O ./SkeletalMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF411NIE/@@download/ENCFF411NIE.bed.gz
wget -O ./FetalSpinalCord_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF252NNW/@@download/ENCFF252NNW.bed.gz
wget -O ./StomachSmoothMuscle_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF777CUK/@@download/ENCFF777CUK.bed.gz
wget -O ./SubcutaneousAbdominalAdiposeTissue_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF467ZOW/@@download/ENCFF467ZOW.bed.gz
wget -O ./TibialNerve_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF221LDE/@@download/ENCFF221LDE.bed.gz
wget -O ./TransverseColon_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF242DHB/@@download/ENCFF242DHB.bed.gz
wget -O ./LeftLungUpperLobe_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF667IAI/@@download/ENCFF667IAI.bed.gz
wget -O ./UrinaryBladder_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF746KLC/@@download/ENCFF746KLC.bed.gz
wget -O ./Uterus_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF045LNJ/@@download/ENCFF045LNJ.bed.gz
wget -O ./Vagina_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF829ZKD/@@download/ENCFF829ZKD.bed.gz
wget -O ./Bcell_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF579EPE/@@download/ENCFF579EPE.bed.gz
wget -O ./Tcell_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF343VOU/@@download/ENCFF343VOU.bed.gz
wget -O ./Monocyte_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF706YSV/@@download/ENCFF706YSV.bed.gz
wget -O ./NKcell_H3K27ac.bed.gz https://www.encodeproject.org/files/ENCFF755QEH/@@download/ENCFF755QEH.bed.gz
while read tissue; do
  echo ${tissue}
  #All peaks
  zcat ${WRKDIR}/data/misc/H3K27ac/${tissue}_H3K27ac.bed.gz | sed 's/^chr//g' | \
  cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/H3K27ac_${tissue}.elements.bed
  #Strong H3K27ac (top 10%)
  cutoff=$( zcat ${WRKDIR}/data/misc/H3K27ac/${tissue}_H3K27ac.bed.gz | cut -f7 | \
  sort -nrk1,1 | perl -e '$d=.1;@l=<>;print $l[int($d*$#l)]' )
  zcat ${WRKDIR}/data/misc/H3K27ac/${tissue}_H3K27ac.bed.gz | awk -v OFS="\t" -v cutoff=${cutoff} \
  '{ if ($7>=cutoff) print $1, $2, $3 }' | sed 's/^chr//g' | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/StrongH3K27ac_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/H3K27ac/*_H3K27ac.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_H3K27ac/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' )
#Conserved H3K27ac
while read tissue; do
  cat ${WRKDIR}/data/master_annotations/noncoding/H3K27ac_${tissue}.elements.bed
done < <( l ${WRKDIR}/data/misc/H3K27ac/*_H3K27ac.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_H3K27ac/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | awk -v OFS="\t" \
  '{ print $1, $2, $3, "0" }' > ${TMPDIR}/H3K27ac_merged.bed
while read tissue; do
  bedtools intersect -c -a ${TMPDIR}/H3K27ac_merged.bed \
  -b ${WRKDIR}/data/master_annotations/noncoding/H3K27ac_${tissue}.elements.bed | \
  awk -v OFS="\t" '{ if ($5>0) $5=1; print $1, $2, $3, $4+$5 }' > \
  ${TMPDIR}/H3K27ac_merged.bed2
  mv ${TMPDIR}/H3K27ac_merged.bed2 ${TMPDIR}/H3K27ac_merged.bed
done < <( l ${WRKDIR}/data/misc/H3K27ac/*_H3K27ac.bed.gz | awk '{ print $9 }' | \
  sed -e 's/_H3K27ac/\t/g' -e 's/\//\t/g' | awk '{ print $(NF-1) }' )
awk -v OFS="\t" '{ if ($4/54>0.5) print $1, $2, $3 }' ${TMPDIR}/H3K27ac_merged.bed > \
${WRKDIR}/data/master_annotations/noncoding/H3K27ac_conserved.elements.bed
awk -v OFS="\t" '{ if ($4/54>0.9) print $1, $2, $3 }' ${TMPDIR}/H3K27ac_merged.bed > \
${WRKDIR}/data/master_annotations/noncoding/H3K27ac_highlyConserved.elements.bed
#TF Binding Sites
cd ${WRKDIR}/data/misc
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
while read TF; do
  echo ${TF}
  #All TF (mild filtering)
  zcat ${WRKDIR}/data/misc/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v OFS="\t" -v TF=${TF} '{ if ($4==TF && $5>200 && $6>1) print $1, $2, $3 }' | \
  sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/TFBS_${TF}.elements.bed
  #Strong TF (heavy filtering)
  zcat ${WRKDIR}/data/misc/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v OFS="\t" -v TF=${TF} '{ if ($4==TF && $5>900 && $6>2) print $1, $2, $3 }' | \
  sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/StrongTFBS_${TF}.elements.bed
done < <( zcat ${WRKDIR}/data/misc/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk '{ if ($6>2) print $4 }' | sort | uniq )
#TFBS union
while read TF; do
  cat ${WRKDIR}/data/master_annotations/noncoding/TFBS_${TF}.elements.bed
done < <( zcat ${WRKDIR}/data/misc/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk '{ if ($6>2) print $4 }' | sort | uniq ) | sort -Vk1,1 -k2,2n -k3,3n | \
  bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/TFBS_union.elements.bed
while read TF; do
  cat ${WRKDIR}/data/master_annotations/noncoding/StrongTFBS_${TF}.elements.bed
done < <( zcat ${WRKDIR}/data/misc/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk '{ if ($6>2) print $4 }' | sort | uniq ) | sort -Vk1,1 -k2,2n -k3,3n | \
  bedtools merge -i - > \
  ${WRKDIR}/data/master_annotations/noncoding/StrongTFBS_union.elements.bed
#Cancer CNV loci (Zack et al, 2013)
##NOTE: manually curated cancer CNV loci from Zack supplement

#Get count of elements per noncoding set (all)
while read list; do
  for dummy in 1; do
    echo -e "${list}"
    #All elements (count)
    cat ${list} | wc -l
    #Write element size to temporary file
    awk '{ print $3-$2 }' ${list} > ${TMPDIR}/element_size.tmp
    #Mean size & std dev
    Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
    cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
    fgrep -v WARNING
    #Autosomal elements (count)
    grep -e '^[0-9]' ${list} | wc -l
    #Write autosomal element size to temporary file
    grep -e '^[0-9]' ${list} | awk '{ print $3-$2 }' > ${TMPDIR}/element_size.tmp
    #Mean autosomal size & std dev
    Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
    options(scipen=1000);\
    cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
    fgrep -v WARNING
  done | paste -s
done < <( l ${WRKDIR}/data/master_annotations/noncoding/*elements.bed | \
  awk '{ print $9 }' | fgrep cancer )
#Get count of elements per noncoding set (ChromHMM, ordered by state)
while read code state; do
  while read list; do
    for dummy in 1; do
      echo -e "${list}"
      #All elements (count)
      cat ${list} | wc -l
      #Write element size to temporary file
      awk '{ print $3-$2 }' ${list} > ${TMPDIR}/element_size.tmp
      #Mean size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
      #Autosomal elements (count)
      grep -e '^[0-9]' ${list} | wc -l
      #Write autosomal element size to temporary file
      grep -e '^[0-9]' ${list} | awk '{ print $3-$2 }' > ${TMPDIR}/element_size.tmp
      #Mean autosomal size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
    done | paste -s
  done < <( l ${WRKDIR}/data/master_annotations/noncoding/*elements.bed | \
    awk '{ print $9 }' | fgrep ${state} )
done < ${WRKDIR}/data/misc/ChromHMM_18way_states.list

#Consolidate all tissue-specific annoation classes into organ system-level tracks
while read tissue; do
  echo ${tissue}
  #CODE:
  # while read anno; do
  #   echo ${anno}
  #   while read file; do
  #     cat ${file}
  #   done < <( awk -v tissue=${tissue} -v anno=${anno} \
  #   '{ if ($1==tissue && $2==anno) print $3 }' \
  #   ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list ) | \
  #   sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - > \
  #   ${WRKDIR}/data/master_annotations/noncoding/${tissue}_MASTER.${anno}.elements.bed
  # done < <( awk -v tissue=${tissue} '{ if ($1==tissue) print $2 }' \
  #   ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
  #   sort | uniq )
  #PARALLELIZE:
  bsub -q short -sla miket_sc -u nobody -J curate_${tissue}_tracks \
  "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_master_organ_group_annotations.sh ${tissue}"
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
sort | uniq )

#Get count of elements per noncoding set (per organ system group)
while read anno; do
  while read tissue; do
    list=${WRKDIR}/data/master_annotations/noncoding/${tissue}_MASTER.${anno}.elements.bed
    for dummy in 1; do
      echo -e "${tissue}\t${anno}"
      echo -e "${list}"
      #All elements (count)
      cat ${list} | wc -l
      #Write element size to temporary file
      awk '{ print $3-$2 }' ${list} > ${TMPDIR}/element_size.tmp
      #Mean size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
      #Autosomal elements (count)
      grep -e '^[0-9]' ${list} | wc -l
      #Write autosomal element size to temporary file
      grep -e '^[0-9]' ${list} | awk '{ print $3-$2 }' > ${TMPDIR}/element_size.tmp
      #Mean autosomal size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      options(scipen=1000);\
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
    done | paste -s
  done < <( awk -v anno=${anno} '{ if ($2==anno) print $1 }' \
    ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
    sort | uniq )
done < <( cut -f2 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
    sort | uniq  )

while read organ; do
  while read list; do
    for dummy in 1; do
      echo ${organ}
      echo -e "${list}"
      #All elements (count)
      cat ${list} | wc -l
      #Write element size to temporary file
      awk '{ print $3-$2 }' ${list} > ${TMPDIR}/element_size.tmp
      #Mean size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
      #Autosomal elements (count)
      grep -e '^[0-9]' ${list} | wc -l
      #Write autosomal element size to temporary file
      grep -e '^[0-9]' ${list} | awk '{ print $3-$2 }' > ${TMPDIR}/element_size.tmp
      #Mean autosomal size & std dev
      Rscript -e "x <- read.table(\"${TMPDIR}/element_size.tmp\",header=F)[,1]; \
      options(scipen=1000);\
      cat(paste(round(mean(x),2),\"\\n\",round(sd(x),2),\"\\n\",sep=\"\"))" | \
      fgrep -v WARNING
    done | paste -s
  done < <( l ${WRKDIR}/data/master_annotations/noncoding/*elements.bed | \
    awk '{ print $9 }' | fgrep ${organ}_MASTER )
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
sort | uniq )

###########################
#####MISC OTHER ANNOTATIONS
###########################
#pcHiC E-P contacts from Javierre et al, Cell, 2016
#Note: must have copied the file ActivePromoterEnhancerLinks.tsv from the paper's supplement
if ! [ -e ${WRKDIR}/data/misc/pcHiC_contacts ]; then
  mkdir ${WRKDIR}/data/misc/pcHiC_contacts
fi
#First: link EP links to gene promoters
sed 's/chr//g' ${WRKDIR}/data/misc/ActivePromoterEnhancerLinks.tsv | sed '1d' | \
bedtools intersect -loj -a - -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.all.bed | \
awk -v OFS="\t" -v FS="\t" '{ print $1, $2, $3, $5, $6, $7, $14, $10 }' | \
awk -v OFS="\t" -v FS="\t" '{ if ($7=="") $7="."; print }' > \
${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.bed
#Second: split EP links based on number of significant
for i in $( seq 1 16 ); do
  # Rscript -e "x <- read.table(\"${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.bed\",header=F);\
  #             out <- x[which(lapply(strsplit(as.character(x[,8]),split=\",\"),length)>=${i}),1:7];\
  #             write.table(out,\"${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_${i}.bed\",\
  #             col.names=F,row.names=F,sep=\"\\t\",quote=F)"
  awk -v OFS="\t" '{ print $4, $5, $6, $7 }' \
  ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_${i}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > \
  ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_${i}.unique_EP_pairs.bed
done
#Affymetrix 6.0 & Illumina Omni SNP array coordinates
cd ${WRKDIR}/data/misc/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snpArrayAffy6.txt.gz
zcat ${WRKDIR}/data/misc/snpArrayAffy6.txt.gz | cut -f2-4 | sed 's/^chr//g' > \
${WRKDIR}/data/misc/Affy6_probes.bed
gzip -f ${WRKDIR}/data/misc/Affy6_probes.bed
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snpArrayIlluminaHumanOmni1_Quad.txt.gz
zcat ${WRKDIR}/data/misc/snpArrayIlluminaHumanOmni1_Quad.txt.gz | cut -f2-4 | sed 's/^chr//g' > \
${WRKDIR}/data/misc/Omni_probes.bed
gzip -f ${WRKDIR}/data/misc/Omni_probes.bed
#Create control probe deserts
while read contig length; do
  paste <( seq 0 10000 $(( ${length}-100000 )) ) \
  <( seq 100000 10000 ${length} ) | \
  awk -v OFS="\t" -v contig=${contig} '{ print contig, $1, $2 }'
done < ${WRKDIR}/data/misc/GRCh37_autosomes.genome > \
${WRKDIR}/data/misc/GRCh37_autosomes.100kbBins_10kbSteps.bed
bedtools intersect -c -a ${WRKDIR}/data/misc/GRCh37_autosomes.100kbBins_10kbSteps.bed \
-b ${WRKDIR}/data/misc/Affy6_probes.bed.gz | awk -v OFS="\t" '{ if ($4<=3) print $1, $2, $3 }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/misc/Affy6_probeDeserts.bed
bedtools intersect -c -a ${WRKDIR}/data/misc/GRCh37_autosomes.100kbBins_10kbSteps.bed \
-b ${WRKDIR}/data/misc/Omni_probes.bed.gz | awk -v OFS="\t" '{ if ($4<=8) print $1, $2, $3 }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/misc/Omni_probeDeserts.bed
cat ${WRKDIR}/data/misc/Affy6_probeDeserts.bed \
${WRKDIR}/data/misc/Omni_probeDeserts.bed | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - > ${WRKDIR}/data/misc/CTRL_probeDeserts.bed
#Disease-associated enhancers (DiseaseEnhancer)
http://biocc.hrbmu.edu.cn/DiseaseEnhancer/RFunctions/diseaseEnh5.1.txt
#Disease-associated enhancers (HEDD)


########################################
#####ANNOTATIONS FOR EXAMPLE LOCUS PLOTS
########################################
if ! [ -e ${WRKDIR}/data/misc/dataForLocusPlots/ ]; then
  mkdir ${WRKDIR}/data/misc/dataForLocusPlots/
fi
#H3K27ac (active enhancers) in adult brain regions. Average across all samples
#81yo male caudate nucleus
#https://www.encodeproject.org/files/ENCFF933PPW/@@download/ENCFF933PPW.tagAlign.gz
#75yo female caudate nucleus
#https://www.encodeproject.org/files/ENCFF855ENP/@@download/ENCFF855ENP.tagAlign.gz
#81yo male cingulate gyrus
#https://www.encodeproject.org/files/ENCFF792CYK/@@download/ENCFF792CYK.tagAlign.gz
#75yo female cingulate gyrus
#https://www.encodeproject.org/files/ENCFF055GDV/@@download/ENCFF055GDV.tagAlign.gz
#81yo male hippocampus
#https://www.encodeproject.org/files/ENCFF085VDH/@@download/ENCFF085VDH.tagAlign.gz
#73yo male hippocampus
#https://www.encodeproject.org/files/ENCFF544ERT/@@download/ENCFF544ERT.tagAlign.gz
#81yo male substantia niagra
#https://www.encodeproject.org/files/ENCFF549OOU/@@download/ENCFF549OOU.tagAlign.gz
for ENCID in ENCFF933PPW ENCFF855ENP ENCFF792CYK ENCFF055GDV ENCFF085VDH ENCFF544ERT ENCFF549OOU; do
  bsub -q filemove -sla miket_sc -J DL_brain_h3k27ac -o /dev/null \
  "cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
  wget https://www.encodeproject.org/files/${ENCID}/@@download/${ENCID}.tagAlign.gz"
done
for ENCID in ENCFF933PPW ENCFF855ENP ENCFF792CYK ENCFF055GDV ENCFF085VDH ENCFF544ERT ENCFF549OOU; do
  zcat ${WRKDIR}/data/misc/dataForLocusPlots/${ENCID}.tagAlign.gz
done | sort -Vk1,1 | bedtools genomecov -bga -i - \
-g <( sed 's/^/chr/g' /data/talkowski/rlc47/src/GRCh37.genome ) > \
${WRKDIR}/data/misc/dataForLocusPlots/adult_brain_H3K27ac_merged.bedGraph
#K3K27me3 (silenced enhancers) in fetal brain (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_h3k27me3 -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF706ZUE/@@download/ENCFF706ZUE.tagAlign.gz;\
ln -fs ./ENCFF706ZUE.tagAlign.gz ./fetal_brain_H3K27me3.tagAlign.gz"
#K3K9me3 (constituative chromatin) in fetal brain (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_h3k9me3 -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF405VIS/@@download/ENCFF405VIS.tagAlign.gz;\
ln -fs ./ENCFF405VIS.tagAlign.gz ./fetal_brain_H3K9me3.tagAlign.gz"
#K3K4me1 (enhancers) in fetal brain (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_h3k4me1 -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF220KAQ/@@download/ENCFF220KAQ.tagAlign.gz;\
ln -fs ./ENCFF220KAQ.tagAlign.gz ./fetal_brain_H3K4me1.tagAlign.gz"
#DHS in fetal brain (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_DHS -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF783LRB/@@download/ENCFF783LRB.bigWig;\
ln -fs ./ENCFF783LRB.bigWig ./fetal_brain_DHS.bigWig"
#DHS in adult frontal cortex (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_DHS -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF420QAT/@@download/ENCFF420QAT.bigWig;\
ln -fs ./ENCFF420QAT.bigWig ./adult_brain_DHS.bigWig"
#HiC heatmap for fetal brain cortical plate (BRAIN)
bsub -q filemove -sla miket_sc -J DL_brain_HiC -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_FBD_IC-heatmap-chr-10k.hdf5.gz;\
ln -fs ./GSE77565_FBD_IC-heatmap-chr-10k.hdf5.gz ./fetal_brain_HiC_10k.hdf5.gz"
# source activate collins_py35
# mkdir ${WRKDIR}/data/misc/dataForLocusPlots/fetalBrainHiC_matrices
# python ${WRKDIR}/bin/utils/hdf2tab/scripts/hdf2tab.py -v \
# --input fetal_brain_HiC_10k.hdf5 \
# -o ${WRKDIR}/data/misc/dataForLocusPlots/fetalBrainHiC_matrices/matrix
# python ${WRKDIR}/bin/utils/hdf2tab/scripts/hdf2tab.py \
# -i fetal_brain_HiC_10k.hdf5 \
# --info
#Fetal brain RNAseq
bsub -q filemove -sla miket_sc -J DL_brain_RNAseq -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF387ESI/@@download/ENCFF387ESI.bigWig;\
ln -fs ./ENCFF387ESI.bigWig ./fetal_brain_RNAseq_minusStrand.bigWig"
bsub -q filemove -sla miket_sc -J DL_brain_RNAseq -o /dev/null \
"cd ${WRKDIR}/data/misc/dataForLocusPlots/; \
wget https://www.encodeproject.org/files/ENCFF862KEW/@@download/ENCFF862KEW.bigWig;\
ln -fs ./ENCFF862KEW.bigWig ./fetal_brain_RNAseq_plusStrand.bigWig"
















