#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for plots of positive control regulatory blocks for per-gene burden tests

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/RegBlockControlSuppFig ]; then
  rm -rf ${WRKDIR}/data/plot_data/RegBlockControlSuppFig
fi
mkdir ${WRKDIR}/data/plot_data/RegBlockControlSuppFig

#####Step 0: Determine eligible regions of the genome to place shuffled segments
#Exclude probe deserts ±50kb
#Exclude gene boundaries from pLI>0.9 genes
#Exclude gene boundaries from dominant disease genes
cat ${WRKDIR}/data/master_annotations/genelists/*dominant* \
    ${WRKDIR}/data/master_annotations/genelists/*Dominant* \
    ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
sort | sed 's/\-/_/g' | fgrep -wf - \
${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -d 50000 -i - | cat - \
<( awk -v OFS="\t" '{ print $1, $2-25000, $3+25000 }' \
   ${WRKDIR}/data/misc/CTRL_probeDeserts.bed | awk -v OFS="\t" \
   '{ if ($2<1) $2=1; print }' ) \
<( grep -e 'X\|Y\|M' ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${TMPDIR}/regBlock_shuffle_blacklist.bed
#Only consider regions where we had at least one qualifying CNV
#Only consider regions with at least one kind of annotation
for pheno in CTRL GERM; do
  zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.DEL.E4.GRCh37.haplosufficient_largeSegmentExcluded_sigGeneExcluded.bed.gz \
       ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.DUP.E4.GRCh37.noncoding_largeSegmentExcluded.bed.gz | \
  fgrep -v "#" | cut -f1-3
done | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${TMPDIR}/regBlock_shuffle_whitelist.tmp
cat ${WRKDIR}/data/master_annotations/noncoding/*elements.bed | \
awk -v OFS="\t" '{ if ($3>$2 && $3-$2>5000) print $1, $2, $3 }' | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools merge -i - | cat - ${TMPDIR}/regBlock_shuffle_whitelist.tmp | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
bedtools subtract -a - -b ${TMPDIR}/regBlock_shuffle_blacklist.bed > \
${TMPDIR}/regBlock_shuffle_whitelist.bed
#Restrict to regions >5kb & remove unplaced contigs
awk '{ if ($3-$2>5000 && $1<=22) print $0 }' ${TMPDIR}/regBlock_shuffle_whitelist.bed | \
fgrep -v random | fgrep -v gl | fgrep -v hap > ${TMPDIR}/regBlock_shuffle_whitelist.bed2
mv ${TMPDIR}/regBlock_shuffle_whitelist.bed2 ${TMPDIR}/regBlock_shuffle_whitelist.bed

#####Step 1a: calcuate distance to nearest constrained gene
#Write lists of protein-coding exons for distance calculations
sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/constrained_exons.bed
#Get empirical distance results
for CNV in DEL DUP; do
  #Oerlap upstream & downstream separately
  paste <( sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | cut -f1-3 | \
           awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' | bedtools closest -D ref -id -io \
           -a - -g /data/talkowski/rlc47/src/GRCh37.genome \
           -b ${TMPDIR}/DDD_exons.bed | awk '{ if ($NF==-1) $NF="NA"; print $NF }' ) \
        <( sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | cut -f1-3 | \
           awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' | bedtools closest -D ref -iu -io \
           -a - -g /data/talkowski/rlc47/src/GRCh37.genome \
           -b ${TMPDIR}/DDD_exons.bed | awk '{ if ($NF==-1) $NF="NA"; print $NF }' ) > \
           ${TMPDIR}/test_distances.txt
  #Quick test: shuffle
  for i in $( seq 1 100 ); do
    #Shuffle
    shuffled=`mktemp`
    bedtools shuffle -g <( grep -e '^[0-9]' /data/talkowski/rlc47/src/GRCh37.genome ) \
    -excl ${TMPDIR}/regBlock_shuffle_blacklist.bed \
    -seed ${i} \
    -i <( sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | cut -f1-3 ) > \
    ${shuffled}
    #Oerlap upstream & downstream separately
    paste <( awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' ${shuffled} | sort -Vk1,1 -k2,2n -k3,3n | \
             bedtools closest -D ref -id -io \
             -a - -g /data/talkowski/rlc47/src/GRCh37.genome \
             -b ${TMPDIR}/DDD_exons.bed | awk '{ if ($NF==-1) $NF="NA"; print $NF }' ) \
          <( awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' ${shuffled} | sort -Vk1,1 -k2,2n -k3,3n | \
             bedtools closest -D ref -iu -io \
             -a - -g <( grep -e '^[0-9]' /data/talkowski/rlc47/src/GRCh37.genome ) \
             -b ${TMPDIR}/DDD_exons.bed | awk '{ if ($NF==-1) $NF="NA"; print $NF }' )
  done > ${TMPDIR}/test_shuf.txt
done 






sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/DDD_exons.bed

sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Dominant_LOF.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/DECIPHER_LOF_exons.bed

sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/DDG2P_HighConf_Dominant_GOF.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/DECIPHER_GOF_exons.bed


sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/BRAIN_MASTER_INTERSECTION.Highly_Expressed.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/BrainExpr_exons.bed


sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_low_confidence.genes.list | \
fgrep -wf - ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | grep -e '^[0-9]' | \
awk -v OFS="\t" '{ print $1, $2, $3, ".", "+" }' > \
${TMPDIR}/ClinGen_HI_exons.bed









