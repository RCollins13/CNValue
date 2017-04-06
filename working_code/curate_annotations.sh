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
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.${filter}.bed
done | paste - - -
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
  #Full path to file
  readlink -f ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.${filter}.bed
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
  awk -v OFS="\t" '{ if ($5>0) print $1, $2, $3, $4+1; else print $1, $2, $3, $4 }' > \
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
#Super Enhancers (dbSUPER)



#PhastCons conservation peaks

#PhyloP conservation peaks

#GERP conservation peaks

#CpG Islands

#Replication timing

#Recombination rate

#Repeat masker

#Get count of elements per noncoding set
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
  awk '{ print $9 }' )
























