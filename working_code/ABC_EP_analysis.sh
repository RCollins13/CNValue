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


#####Intersect cleaned enhancers with cleaned deletion regulatory blocks
CNV=DEL; VF=E4
#Get list of high-confidence deletion regulatory blocks
while read ID sites; do
  for wrapper in 1; do
    echo -e "${ID}\t${CNV}"
    #Any inheritance
    for mem in probands siblings; do
      echo -e "${sites}" | sed -e 's/\_/\t/g' -e 's/\;/\n/g' | \
      bedtools intersect -u -b - \
      -a ${TMPDIR}/SSC_CNVs.p10E9.hg19.${mem}.${CNV}.bed | \
      cut -f4 | sort | uniq | wc -l
    done
    #De novo CNVs only
    for mem in probands siblings; do
      echo -e "${sites}" | sed -e 's/\_/\t/g' -e 's/\;/\n/g' | \
      bedtools intersect -u -b - \
      -a ${TMPDIR}/SSC_CNVs.p10E9.hg19.${mem}.${CNV}.bed | fgrep DeNovo | \
      cut -f4 | sort | uniq | wc -l
    done
  done | paste -s
done < <( sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
          cut -f4,5 ) | awk -v OFS="\t" \
'{ if ($3>=$4 && $5>=$6 && $3<10) print $1 }' | \
fgrep -wf - ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
awk -v OFS="\t" '{ if ($6<10 && $7>5) print $0 }' | cut -f1-4 | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/HC_DEL_regBlocks.bed
#Print list of hotspots with enhancers of constrained or disease-associated genes
bedtools intersect -wa -wb -a ${TMPDIR}/HC_DEL_regBlocks.bed \
-b ${WRKDIR}/data/misc/ABC_EP/RLC_processed/ABC_EP.merged_tissues.10kb_max.bed | \
fgrep -wf <( cat ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list \
                 ${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Dominant_LOF.genes.list \
                 ${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_low_confidence.genes.list | \
                 sort | uniq )


#####Estimate expected number of hits by chance
#Build list of excluded loci
cat <( awk -v OFS="\t" '{ print $1, $2-40000, $3+40000 }' \
       ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
       sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ if ($2<0) $2=0; print $0 }') \
${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > ${TMPDIR}/excluded_loci.tmp.bed
#Subset enhancers to genes of interest and not within 40kb of an exon of a non-haplosufficient gene
cat ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list \
${WRKDIR}/data/master_annotations/genelists/DDG2P_AnyConf_Dominant_LOF.genes.list \
${WRKDIR}/data/master_annotations/genelists/ClinGen_haploinsufficient_low_confidence.genes.list | \
sort | uniq | fgrep -wf - \
${WRKDIR}/data/misc/ABC_EP/RLC_processed/ABC_EP.merged_tissues.10kb_max.bed | \
bedtools intersect -v -a - -b ${TMPDIR}/excluded_loci.tmp.bed | cut -f1-3 > \
${TMPDIR}/eligible_enhancers.bed
#Get observed overlap
bedtools intersect -u \
-a ${TMPDIR}/HC_DEL_regBlocks.bed \
-b ${TMPDIR}/eligible_enhancers.bed | wc -l
#Shuffle 100 times to estimate expected overlap
for i in $( seq 1 100 ); do
  bedtools shuffle -noOverlapping -seed ${i} \
  -excl ${TMPDIR}/excluded_loci.tmp.bed \
  -i ${TMPDIR}/HC_DEL_regBlocks.bed \
  -g /data/talkowski/rlc47/src/GRCh37.genome | \
  bedtools intersect -u -a - \
  -b ${TMPDIR}/eligible_enhancers.bed | wc -l
done






