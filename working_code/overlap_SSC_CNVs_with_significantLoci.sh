#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to find SSC probands with CNVs overlapping significant genes & reg blocks

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Make proband files for hg19 DEL & DUP with sample IDs
#Subset CNVs with pCNV ≤ 1x10^-9; print chr, start, end, ID, CN, CNV, inheritance, pCNV
awk -v OFS="\t" '{ if ($4=="Proband" && $45<=0.000000001) print $9, $10, $11, $7, $29, $43, $44, $45 }' \
/data/talkowski/Samples/SFARI/OLD/CNV_comparison/SSC_SNP_Genotyping_CNVs.txt | \
fgrep -v "No_CNV" | awk '{ if ($5!=2 && $5!="2?") print $0 }' > ${TMPDIR}/SSC_CNVs.p10E9.hg18.probands.bed
#Lift over to hg19
liftOver -minMatch=0.5 -bedPlus=5 ${TMPDIR}/SSC_CNVs.p10E9.hg18.probands.bed \
/data/talkowski/rlc47/src/hg18ToHg19.over.chain ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed \
${TMPDIR}/SSC_CNVs.p10E9.hg18_liftFail.bed
sed -e 's/chr//g' -e 's/\?//g' ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed2
mv ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed2 ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed
#Split by del/dup
awk -v OFS="\t" '{ if ($5<2) print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed > \
${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.DEL.bed
awk -v OFS="\t" '{ if ($5>2) print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.bed > \
${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.DUP.bed

#####Make master list of all SSC overlapping CNVs
for wrapper in 1; do
  #Intersect vs. genes
  VF=E4; context=exonic
  for CNV in DEL DUP; do
    fgrep -wf \
    <( cut -f1 ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | sed 's/\-/_/g' ) \
    <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed ) | \
    bedtools intersect -wb -a - -b ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.${CNV}.bed | \
    awk -v OFS="\t" '{ print $5, $6, $7, $4, $8, $10, $11 }'
  done
  #Intersect vs. regulatory blocks
  for CNV in DEL DUP; do
    bedtools intersect -wb \
    -a <( cut -f1-4 ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | sed '1d' ) \
    -b ${TMPDIR}/SSC_CNVs.p10E9.hg19.probands.${CNV}.bed | \
    awk -v OFS="\t" '{ print $5, $6, $7, $4, $8, $10, $11 }'
  done
done | sort -Vk1,1 -Vk4,4 -k6,6 -k2,2n -k3,3n -Vk5,5 -k7,7 | uniq | \
awk -v OFS="\t" '{ if ($3-$2>50000 && $3-$2<5000000) print $0 }' > \
${TMPDIR}/SSC_proband_CNVs_overlapping_sigLoci.bed