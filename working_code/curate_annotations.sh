#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master annotation curation code
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
#Haplosufficient genes (pLI < 0.1)
awk '{ if ($20<0.1 && $20!="NA") print $2 }' \
${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | sort | uniq > \
${WRKDIR}/data/master_annotations/genelists/ExAC_haplosufficient.genes.list
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
  awk '{ print $9 }' ) | paste - - -



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
<( sed '1d' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_CER.bed | cut -f1,4,5 ) | \
sort -Vk1,1 -k2,2n -k3,3n | awk '{ if ($3-$2<=5000000) print $0 }' | bedtools merge -i - | \
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
done < <( tail -n3 ${WRKDIR}/data/misc/SuperEnhancer_tissues.list )
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
  awk -v OFS="\t" '{ if ($4>=1) print $1, $2, $3 }' \
  ${WRKDIR}/data/misc/PhyloP/chr${chr}.phyloP46way.placental.bg | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 100 -i - | \
  awk -v OFS="\t" '{ if ($3-$2>=200) print $1, $2, $3 }' | sed 's/^chr//g'
done > \
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
#Recombination rate

#Repeat masker

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
sed 's/\ /_/g' ${WRKDIR}/data/misc/gwas_catalog_v1.0.1-associations_e88_r2017-04-03.tsv | \
awk -v FS="\t" -v OFS="\t" '{ if ($34=="N") print $0 }' | awk -v FS="\t" -v OFS="\t" \
'$1 ~ /2013|2014|2015|2016|2017/ { print $12, $13, $13+1, $35 }' | \
awk -v OFS="\t" -v FS="\t" '{ if ($1!="" && $2!="" && $3!="" && $4!="") print $0 }' | cut -f4 | sort | uniq

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
  awk '{ print $9 }' | fgrep PMD )
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























