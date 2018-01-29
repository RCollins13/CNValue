#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to process & analyze all ABC EP connections from Lander lab

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Prep directories
for dir in ABC_EP ABC_EP/RLC_processed; do
  if ! [ -e ${WRKDIR}/data/misc/${dir} ]; then
    mkdir ${WRKDIR}/data/misc/${dir}
  fi
done

#####Rsync data
bsub -q filemove -sla miket_sc -J ABC_EP_rsync \
"rsync -ravz \
rcollins@gold.broadinstitute.org:/seq/lincRNA/RAP/ABC/171116_pred/* \
${WRKDIR}/data/misc/ABC_EP/"

#####Process all tracks per tissue
#Write list of all tissues
ls -ltrh ${WRKDIR}/data/misc/ABC_EP/ | awk '{ print $9 }' | sed '1d' | \
grep -wv 'analysis\|qc\|RLC_processed' > \
${WRKDIR}/data/misc/ABC_EP/tissues.list
#Write simplified list of distal enhancers to file per tissue
while read tissue; do
  if [ -e ${WRKDIR}/data/misc/ABC_EP/${tissue}/EnhancerPredictions.bed ]; then
    echo ${tissue}
    sed '1d' ${WRKDIR}/data/misc/ABC_EP/${tissue}/EnhancerPredictions.bed | \
    cut -f4 | fgrep "distal" | sed 's/|/\t/g' | cut -f2 | \
    sed -e 's/:\|>\|\-/\t/g' | sed 's/[\t]+/\t/g' | sed 's/^chr//g' | \
    awk -v OFS="\t" -v tissue=${tissue} '{ print $1, $2, $3, $4, tissue }' > \
    ${WRKDIR}/data/misc/ABC_EP/RLC_processed/${tissue}.enhancers.bed
  fi
done < ${WRKDIR}/data/misc/ABC_EP/tissues.list
#Merge distal enhancers across all tissues into master file
while read tissue; do
  if [ -e ${WRKDIR}/data/misc/ABC_EP/${tissue}/EnhancerPredictions.bed ]; then
    awk -v FS="\t" '{ if ($3-$2<10000) print $0 }' \
    ${WRKDIR}/data/misc/ABC_EP/RLC_processed/${tissue}.enhancers.bed
  fi
done < ${WRKDIR}/data/misc/ABC_EP/tissues.list | \
sort -Vk1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
bedtools merge -c 4,5 -o distinct > \
${WRKDIR}/data/misc/ABC_EP/RLC_processed/ABC_EP.merged_tissues.10kb_max.bed


#####Intersect cleaned enhancers with significant regulatory blocks
#Print genes with enhancers 
for CNV in DEL DUP; do
  sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
  cut -f1-4 | bedtools intersect -wa -wb -a - \
  -b ${WRKDIR}/data/misc/ABC_EP/RLC_processed/ABC_EP.merged_tissues.10kb_max.bed
done











