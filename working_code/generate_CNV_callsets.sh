#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to create all master CNV sets for rCNVmap project

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Gather TCGA CNVs
#Determine number of unique tumor types represented (samples, not CNVs)
awk '{ print "s/"$1"/"$3"/g" }' ${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.txt > \
${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.sed
sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | cut -f1 | sort | uniq | \
awk -v FS="-" '{ print $2 }' | sed -f ${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.sed | \
sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }'
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs
#Filter CNV classes by log2 ratios
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CTRL.${CNV}.raw.bed
done
#Cancers
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed
done
paste <( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF<=-1) print $0 }' | \
 awk -v OFS="\t" '{ print $2, $3, $4, "TCGA_CNCR_DEL", "DEL" }' ) \
<( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF<=-1) print $1 }' | \
awk -v FS="-" '{ print $2 }' | sed -f ${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.sed | \
awk '{ print toupper($0) }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, "24071852" }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.DEL.raw.bed
paste <( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF>=0.5849625) print $0 }' | \
 awk -v OFS="\t" '{ print $2, $3, $4, "TCGA_CNCR_DUP", "DUP" }' ) \
<( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF>=0.5849625) print $1 }' | \
awk -v FS="-" '{ print $2 }' | sed -f ${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.sed | \
awk '{ print toupper($0) }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, "24071852" }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.DUP.raw.bed
#Clean controls
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CTRL.${CNV}.raw.bed
done
cat ${TMPDIR}/../misc_CNVs/TCGA_normals/*txt | fgrep -v Chromosome | \
awk -v OFS="\t" '{ if ($NF<=-1) print $2, $3, $4, "TCGA_CTRL_DEL_", "DEL", "HEALTHY_CONTROL", "24071852" }' | \
sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ print $1, $2, $3, $4""NR, $5, $6, $7 }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CTRL.DEL.raw.bed
cat ${TMPDIR}/../misc_CNVs/TCGA_normals/*txt | fgrep -v Chromosome | \
awk -v OFS="\t" '{ if ($NF>=0.5849625) print $2, $3, $4, "TCGA_CTRL_DUP_", "DUP", "HEALTHY_CONTROL", "24071852" }' | \
sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ print $1, $2, $3, $4""NR, $5, $6, $7 }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CTRL.DUP.raw.bed
gzip -f ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/*.raw.bed

#####Gather raw SSC CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs
#Subset CNVs with pCNV ≤ 1x10^-9; print chr, start, end, ID, CN, CNV, inheritance, pCNV
awk -v OFS="\t" '{ if ($45<=0.000000001) print $9, $10, $11, $7, $29, $43, $44, $45 }' \
/data/talkowski/Samples/SFARI/OLD/CNV_comparison/SSC_SNP_Genotyping_CNVs.txt | \
fgrep -v "No_CNV" | awk '{ if ($5!=2 && $5!="2?") print $0 }' > ${TMPDIR}/SSC_CNVs.p10E9.hg18.bed
#Lift over to hg19
liftOver -minMatch=0.5 -bedPlus=5 ${TMPDIR}/SSC_CNVs.p10E9.hg18.bed \
/data/talkowski/rlc47/src/hg18ToHg19.over.chain ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed \
${TMPDIR}/SSC_CNVs.p10E9.hg18_liftFail.bed
sed 's/chr//g' ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed > ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed2
mv ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed2 ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed
#Split by del/dup
awk -v OFS="\t" '{ if ($5<2) print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed > \
${TMPDIR}/SSC_CNVs.p10E9.hg19.DEL.bed
awk -v OFS="\t" '{ if ($5>2) print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.bed > \
${TMPDIR}/SSC_CNVs.p10E9.hg19.DUP.bed
#Make ID matching file for phenotypes
fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/SSC_phenotype_key.txt | \
awk '{ print "s/"$1"/"$2"/g" }' > ${WRKDIR}/bin/rCNVmap/misc/SSC_phenotype_key.sed
#Split by case/control
for CNV in DEL DUP; do
  #Controls
  cat <( echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" ) \
  <( awk '$4 ~ /fa|mo/ { print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | sed 's/\?//g' | \
  awk -v CNV=${CNV} -v OFS="\t" '{ print $1, $2, $3, "SSC_CTRL_"CNV"_"NR, CNV, "HEALTHY_CONTROL", "26402605" }' ) > \
  ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_CTRL.${CNV}.raw.bed
  #Cases
  cat <( echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" ) > \
  ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_GERM.${CNV}.raw.bed
  paste <( awk '$4 ~ /p/ { print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.${CNV}.bed | \
  awk -v CNV=${CNV} -v OFS="\t" '{ print $1, $2, $3, "SSC_GERM_"CNV"_", CNV }' ) \
  <( awk '$4 ~ /p/ { print $4 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.${CNV}.bed | \
  sed -f ${WRKDIR}/bin/rCNVmap/misc/SSC_phenotype_key.sed ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ print $1, $2, $3, $4""NR, $5, $6, "26402605" }' |  \
  sed 's/\?//g' >> ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_GERM.${CNV}.raw.bed
done
gzip -f ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/*.bed

#####Gather all raw DGV CNVs (including studies not used in analysis)
while read study; do
  echo "${study}"
  fgrep -w "${study}" /scratch/miket/rlc47temp/DGV_CNVs_all/GRCh37_hg19_variants_2016-05-15.txt > \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}.calls.txt
done < <( sed '1d' /scratch/miket/rlc47temp/DGV_CNVs_all/GRCh37_hg19_variants_2016-05-15.txt | cut -f7 | sort | uniq )
while read study; do
  echo ${study}
  rm -rf /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs
  if ! [ -e /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs ]; then
    mkdir /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs
  fi
  nstudy=$( echo ${study} | sed 's/_et_al/\t/g' | cut -f1 )
  #Deletions
  echo -e "chr\tstart\tend\tVarID\tCNV\tObservations\tSampSize\tFreq\tPMID" > \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs/${nstudy}.DEL.raw.bed
  awk -v FS="\t" -v OFS="\t" '{ if ($17>0 && ($6=="deletion" || $6=="loss")) print $2, $3, $4, $1, "DEL", $17, $15, $14, $8 }' \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}.calls.txt | sort -Vk1,1 -k2,2n -k3,3n >> \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs/${nstudy}.DEL.raw.bed
  #Duplications
  echo -e "chr\tstart\tend\tVarID\tCNV\tObservations\tSampSize\tFreq\tPMID" > \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs/${nstudy}.DUP.raw.bed
  awk -v FS="\t" -v OFS="\t" '{ if ($16>0 && ($6=="duplication" || $6=="gain" || $6=="insertion" || $6=="tandem duplication")) print $2, $3, $4, $1, "DUP", $16, $15, $14, $8 }' \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}.calls.txt | sort -Vk1,1 -k2,2n -k3,3n > \
  /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs/${nstudy}.DUP.raw.bed  
  #Gzip
  gzip /scratch/miket/rlc47temp/DGV_CNVs_all/${study}_CNVs/*bed*
done < <( sed '1d' /scratch/miket/rlc47temp/DGV_CNVs_all/GRCh37_hg19_variants_2016-05-15.txt | \
  cut -f7 | sort | uniq )

#NOTE: ITSARA REPLACED BY COOPER ET AL 2011
# #####Gather Itsara CNVs
# if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs ]; then
#   rm -r ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs
# fi
# mkdir ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs
# for CNV in DEL DUP; do
#   echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
#   ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/Itsara_CTRL.${CNV}.raw.bed
#   while read chr start end sID skip n total PMID; do
#     for times in $( seq 1 ${n} ); do
#       echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tHEALTHY_CONTROL\t${PMID}"
#     done
#   done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Itsara_CNVs/Itsara.${CNV}.raw.bed.gz | sed '1d' ) | \
#   sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
#   '{ print $1, $2, $3, "Itsara_CTRL_"CNV"_"NR, $5, $6 }' >> \
#   ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/Itsara_CTRL.${CNV}.raw.bed
# done
# gzip -f ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/*.raw.bed

####Gather Cooper control CNVs
#Clean Cooper CNVs
sed 's/\",/\t/g' ${TMPDIR}/../misc_CNVs/supporting_variants_for_nstd54.csv | fgrep hg19 | \
sed '1d' | tr -d "\"" | awk -v OFS="\t" -v FS="\t" \
'{ if ($5>1) print $14, $17, $18, "ID", $4, "HEALTHY_CONTROL", "21841781" }' | \
sed 's/copy\ number\ //g' | sed -e 's/gain/DUP/g' -e 's/loss/DEL/g' > \
${TMPDIR}/../misc_CNVs/Cooper_2011_control_CNVcalls.GRCh37.bed
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs/Cooper_CTRL.${CNV}.raw.bed
  fgrep -w ${CNV} ${TMPDIR}/../misc_CNVs/Cooper_2011_control_CNVcalls.GRCh37.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Cooper_CTRL_"CNV"_"NR, $5, $6, $7 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs/Cooper_CTRL.${CNV}.raw.bed
done
gzip -f ${WRKDIR}/data/CNV/CNV_RAW/Cooper_CNVs/Cooper_CTRL.*.bed

#####Gather Shaikh CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs/Shaikh_CTRL.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tHEALTHY_CONTROL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Shaikh_CNVs/Shaikh.${CNV}.raw.bed.gz | sed '1d' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Shaikh_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs/Shaikh_CTRL.${CNV}.raw.bed
done
gzip -f ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs/*.raw.bed

#####Gather Suktitipat CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs/Suktitipat_CTRL.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tHEALTHY_CONTROL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Suktitipat_Thai_CNVs/Suktitipat.${CNV}.raw.bed.gz | sed '1d' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Suktitipat_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs/Suktitipat_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Suktitipat_CNVs/*.raw.bed

#####Gather Uddin CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs/Uddin_CTRL.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tHEALTHY_CONTROL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Uddin_Ontario_CNVs/Uddin.${CNV}.raw.bed.gz | sed '1d' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Uddin_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs/Uddin_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Uddin_CNVs/*.raw.bed

#####Gather Vogler CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs/Vogler_CTRL.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tHEALTHY_CONTROL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Vogler_CNVs/Vogler.${CNV}.raw.bed.gz | sed '1d' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Vogler_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs/Vogler_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs/*.raw.bed

# #####Gather Huang CNVs
# if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs ]; then
#   rm -r ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs
# fi
# mkdir ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs
# #Headers
# for CNV in DEL DUP; do
#   echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
#   ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs/Huang_GERM.${CNV}.raw.bed
# done
# #DELs
# awk -v OFS="\t" '{ if ($3==2 && $10>2) print $6, $7, $8, "Huang_GERM_DEL", "DEL", "TOURETTES", "28641109" }' \
# ${TMPDIR}/Huang_2017_CNVs.raw.bed | sort -Vk1,1 -k2,2n -k3,3n | \
# awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, $7 }' >> \
# ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs/Huang_GERM.DEL.raw.bed
# #DUPs
# awk -v OFS="\t" '{ if ($3==2 && $10<2) print $6, $7, $8, "Huang_GERM_DUP", "DUP", "TOURETTES", "28641109" }' \
# ${TMPDIR}/Huang_2017_CNVs.raw.bed | sort -Vk1,1 -k2,2n -k3,3n | \
# awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, $7 }' >> \
# ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs/Huang_GERM.DUP.raw.bed
# gzip ${WRKDIR}/data/CNV/CNV_RAW/Huang_CNVs/*.raw.bed

#####Gather Talkowski SigGen CNVs
#Download from database (MySQL code below)
# SELECT `Chr`, `Start`, `Stop`, `Diag_state`, `ID`, `Source`, `Indication`, `Ref` \
# FROM all_cnvs WHERE Source = 'SIGNATUREGENOMICS' OR Source = 'SIGNATUREGENOMICS2' OR Source = 'SIGNATUREGENOMICS3';
#Clean file
sed -e 's/\ /_/g' -e 's/,/\;/g' ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.bed2
mv ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.bed2 ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.bed
#Split by native reference
for ref in hg18 hg19; do
  awk -v ref=${ref} '{ if ($NF==ref) print $0 }' ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.bed > \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.${ref}_native.bed
done
#Liftover hg18 native & merge with hg19 native
liftOver -minMatch=0.5 -bedPlus=5 \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg18_native.bed \
/data/talkowski/rlc47/src/hg18ToHg19.over.chain \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_lifted.bed \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg18_liftFail.bed
cat ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_lifted.bed \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_native.bed | sed 's/^chr//g' > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed
#Benchmarking: titrate number of unique cases in our database with a CNV with one and only one hit in Eichler's sets
for r in $( seq 0.05 0.05 1 ); do
  echo ${r}
  bedtools intersect -r -f ${r} -c \
  -a ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed \
  -b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
  awk '{ if ($NF==1) print $5 }' | sort | uniq | wc -l
done | paste - -
#Get pairs of subject IDs and phenotypes of cases that have at least one perfect match
bedtools intersect -r -f 1.0 -c \
-a ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed \
-b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
awk -v OFS="\t" '{ if ($NF==1) print $1, $2, $3, $4, $5, $7 }' | \
bedtools intersect -r -f 1.0 -wa -wb -a - \
-b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
awk -v OFS="\t" '{ print $5, $12, $6, $13 }' | sort -nk1,1 | uniq > \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt
cut -f1 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sort | uniq -c | \
awk -v OFS="\t" '{ if ($1==1) print $2 }' | fgrep -wf - \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt > \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2
mv ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt
cut -f2 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sort | uniq -c | \
awk -v OFS="\t" '{ if ($1==1) print $2 }' | fgrep -wf - \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt > \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2
mv ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt
sed -e 's/\_\;/\;/g' -e 's/\;\_/\;/g' ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt > \
${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2
mv ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt2 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt
#Summarize linked phenotypes - Eichler
cut -f4 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sort | uniq -c | \
awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1
cut -f4 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sed 's/\;/\n/g' | sort | uniq -c | \
awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1
#Summarize linked phenotypes - Talkowski
cut -f3 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sort | uniq -c | \
awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1
# cut -f3 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sort | uniq -c | \
# awk -v OFS="\t" '{ if ($1<10) print $1 }' | awk '{ sum+=$1 }END{ print sum }'
cut -f3 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sed 's/\;/\n/g' | \
sed 's/^_//g' | sort | uniq -c | awk -v OFS="\t" '{ if ($1>9) print $1, $2 }' | sort -nrk1,1
cut -f3 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | sed 's/\;/\n/g' | \
sed 's/^_//g' | sort | uniq -c | awk -v OFS="\t" '{ if ($1<=9) print $2 }' | wc -l

#Benchmarking: titrate number of unique cases in our database with no CNVs with at least one hit in Eichler's sets
for r in $( seq 0.05 0.05 1 ); do
  echo ${r}
  bedtools intersect -r -f ${r} -c \
  -a ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed \
  -b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
  awk '{ if ($NF==1) print $5 }' | sort | uniq | \
  fgrep -wvf - <( cut -f5 ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed ) | \
  sort | uniq | wc -l
done | paste - -
bedtools intersect -r -f ${r} -c \
  -a ${TMPDIR}/../misc_CNVs/TGDB_CNVs.SigGen_ALL.hg19_merged.bed \
  -b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
  awk '{ if ($NF==1) print $0 }' | \
  bedtools intersect -r -f ${r} -wa -wb -a - \
  -b ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | less -S

# #Split by CNV class
# sed 's/\ /_/g' ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed2
# mv ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed2 \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed
# awk -v OFS="\t" -v FS="\t" '{ if ($4=="0,200,0" || $4=="COPY_LOSS" || \
#   $4=="COPY_LOSS." || $4=="COPY_LOSSETION" || $4=="COPY_LOSS_(MOS)" || \
#   $4=="kiss" || $4=="los" || $4=="loss" || $4=="loss_" || $4=="Loss" || \
#   $4=="LOSS" || $4=="MOS_COPY_LOSS" || $4=="NULL_COPY_LOSS") print $0 }' \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.DEL.bed
# awk -v OFS="\t" -v FS="\t" '{ if ($4=="220,0,0" || $4=="255,0,0" || $4=="AMP" || \
#   $4=="AMPLIFICATION" || $4=="COPY_GAIN" || $4=="COPY_GAIN_" || \
#   $4=="COPY_GAIN_(MOS?)" || $4=="COPY_GAIN_(MOS)" || $4=="dup" || $4=="g" || \
#   $4=="gain" || $4=="gain_" || $4=="Gain" || $4=="GAIN" || \
#   $4=="MOSAIC_COPY_GAIN" || $4=="MOS_COPY_GAIN") print $0 }' \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
# ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.DUP.bed
# #Split by native reference
# for CNV in DEL DUP; do
#   for ref in hg18 hg19; do
#     awk -v ref=${ref} '{ if ($NF==ref) print $0 }' \
#     ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.bed > \
#     ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.${ref}_native.bed
#   done
# done
# #Liftover to hg19
# for CNV in DEL DUP; do
#   liftOver -minMatch=0.5 -bedPlus=4 \
#   ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg18_native.bed \
#   /data/talkowski/rlc47/src/hg18ToHg19.over.chain \
#   ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_lifted.bed \
#   ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg18_liftFail.bed
#   cat ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_native.bed \
#   ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_lifted.bed | \
#   sed 's/chr//g' | sort -Vk1,1 -k2,2n -k3,3n > \
#   ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_merged.bed
# done

#####Gather Coe CNVs
#Clean Coe control & case data
sed 's/\",/\t/g' ${TMPDIR}/../misc_CNVs/supporting_variants_for_nstd100.csv | fgrep hg19 | \
sed '1d' | tr -d "\"" | awk -v OFS="\t" -v FS="\t" \
'{ if ($5>1) print $14, $17, $18, "ID", $4, "HEALTHY_CONTROL", "25217958" }' | \
sed 's/copy\ number\ //g' | sed -e 's/gain/DUP/g' -e 's/loss/DEL/g' > \
${TMPDIR}/../misc_CNVs/Coe_2014_control_CNVcalls.GRCh37.bed
sed 's/\",/\t/g' ${TMPDIR}/../misc_CNVs/supporting_variants_for_nstd100.csv | fgrep hg19 | \
sed '1d' | tr -d "\"" | awk -v OFS="\t" -v FS="\t" \
'{ if ($6=="Oligo aCGH") print $14, $17, $18, "ID", $4, $9, "25217958" }' | \
sed 's/copy\ number\ //g' | sed -e 's/gain/DUP/g' -e 's/loss/DEL/g' > \
${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.GRCh37.bed
#Curate Coe & Cooper phenotype data
#NOTE: since Cooper has more detailed phenotypes, merge the two and take nonredunant phenos
fgrep -v "Study ID" ${TMPDIR}/../misc_CNVs/samples_for_nstd54_1.csv | \
sed 's/\ /_/g' | sed 's/\",/\t/g' | tr -d "\"" | \
awk -v FS="\t" -v OFS="\t" '{ print $3, toupper($9) }' | sed 's/,_/\;/g' > \
${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list
fgrep -wvf <( cut -f1 ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list ) \
${TMPDIR}/../misc_CNVs/samples_for_nstd100_1.csv | fgrep -v "Study ID" | \
sed 's/\ /_/g' | sed 's/\",/\t/g' | tr -d "\"" | awk -v FS="\t" -v OFS="\t" \
'{ print $3, toupper($9) }' | sed 's/,_/\;/g' >> \
${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list
#Replace phenotypes for cases able to be uniquely identified from overlaps with our CNVs
cat <( cut -f2 ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | \
  fgrep -wvf - ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list ) \
<( awk -v OFS="\t" '{ print $2, $3";"$4 }' ${TMPDIR}/SigGen_UniqueMatches.100pct_recip.txt | \
  sed -e 's/\-//g' -e 's/\?//g' -e 's/\///g') | \
sort -nk1,1 | awk -v OFS="\t" '{ print $1, $2 }' > ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list2
mv ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list2 \
${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list
#Prepare sed file
awk '{ print "s/"$1"/"$2"/g" }' ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list > \
${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.sed
#Assign phenotypes to Coe cases
paste <( cut -f1-6 ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.GRCh37.bed ) \
<( cut -f6 ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.GRCh37.bed | \
sed -f ${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.sed ) | \
awk -v OFS="\t" '{ print $0, "25217958" }' > \
${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed
#Curate Coe raw CNV files
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs
for CNV in DEL DUP; do
  for pheno in GERM CTRL; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe_${pheno}.${CNV}.raw.bed
  done
  fgrep -w ${CNV} ${TMPDIR}/../misc_CNVs/Coe_2014_case_CNVcalls.wPhenos.GRCh37.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Coe_GERM_"CNV"_"NR, $5, $7, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe_GERM.${CNV}.raw.bed
  fgrep -w ${CNV} ${TMPDIR}/../misc_CNVs/Coe_2014_control_CNVcalls.GRCh37.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Coe_CTRL_"CNV"_"NR, $5, $6, $7 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe_CTRL.${CNV}.raw.bed
done
gzip -f ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe*.bed

#####Gather PGC CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs
#Liftover to hg19
liftOver -minMatch=0.5 -bedPlus=3 ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.txt \
/data/talkowski/rlc47/src/hg18ToHg19.over.chain \
${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed \
${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg18_liftFail.bed
sed 's/chr//g' ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed > \
${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed2
mv ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed2 \
${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed
#Split by CNV class
fgrep -v "#" ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed | fgrep del | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.DEL.bed
fgrep -v "#" ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.bed | fgrep dup | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.DUP.bed
#Split by case/control
for CNV in DEL DUP; do
  #Case (SCZ)
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_GERM.${CNV}.raw.bed
  fgrep -w Case ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "PGC_GERM_"CNV"_"NR, CNV, "SCHIZOPHRENIA", "27869829" }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_GERM.${CNV}.raw.bed
  #Control
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_CTRL.${CNV}.raw.bed
  fgrep -w Control ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "PGC_CTRL_"CNV"_"NR, CNV, "HEALTHY_CONTROL", "27869829" }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/*.raw.bed

#####Gather Talkowski CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
#Download CNVs from CHB and GeneDX from database (MySQL code below)
# SELECT `Chr`, `Start`, `Stop`, `Diag_state`, `ID`, `Source`, `Indication`, `Ref` \
# FROM all_cnvs WHERE Source = 'CHB' OR Source = 'GENEDX';
#Clean phenotype text
sed '1d' ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed | \
sed -e 's/COORDINATES\ ALTERED\ CONTACT\ CHB//g' -e 's/4\ COPIES//g' | \
sed -e 's/\ /_/g' -e 's/,/\;/g' -e 's/_\;/\;/g' -e 's/\;\_/\;/g' -e 's/(/_/g' \
  -e 's/)/_/g'  -e 's/_\t/\t/g' -e 's/\.//g' -e 's/\-//g' -e 's/\?//g' -e 's/\"//g' | \
sed -e 's/\;\t/\t/g' -e 's/\;\;\t/\t/g' -e 's/_\t/\t/g' | sed -e 's/_\t/\t/g' | \
awk -v OFS="\t" -v FS="\t" '{ print $1, $2, $3, $4, $6"_"$5, $7, $8 }' | sed -e 's/\;\t/\t/g' -e 's/\t__//g' | \
awk -v OFS="\t" -v FS="\t" '{ if ($6=="") $6="NOT_INDICATED"; print }' > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed2
mv ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed2 ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed
#Extract unique ID+phenotype pairings
cut -f5,6 ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed | sort | uniq > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.phenotypes.list
#Split by CNV class
awk -v OFS="\t" -v FS="\t" '{ if ($4=="0,200,0" || $4=="COPY_LOSS" || \
  $4=="COPY_LOSS." || $4=="COPY_LOSSETION" || $4=="COPY_LOSS_(MOS)" || \
  $4=="kiss" || $4=="los" || $4=="loss" || $4=="loss_" || $4=="Loss" || \
  $4=="LOSS" || $4=="MOS_COPY_LOSS" || $4=="NULL_COPY_LOSS") print $0 }' \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.DEL.bed
awk -v OFS="\t" -v FS="\t" '{ if ($4=="220,0,0" || $4=="255,0,0" || $4=="AMP" || \
  $4=="AMPLIFICATION" || $4=="COPY_GAIN" || $4=="COPY_GAIN_" || \
  $4=="COPY_GAIN_(MOS?)" || $4=="COPY_GAIN_(MOS)" || $4=="dup" || $4=="g" || \
  $4=="gain" || $4=="gain_" || $4=="Gain" || $4=="GAIN" || \
  $4=="MOSAIC_COPY_GAIN" || $4=="MOS_COPY_GAIN") print $0 }' \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.DUP.bed
#Split by native reference
for CNV in DEL DUP; do
  for ref in hg18 hg19; do
    awk -v ref=${ref} '{ if ($NF==ref) print $0 }' \
    ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.bed > \
    ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.${ref}_native.bed
  done
done
#Liftover to hg19
for CNV in DEL DUP; do
  liftOver -minMatch=0.5 -bedPlus=4 \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg18_native.bed \
  /data/talkowski/rlc47/src/hg18ToHg19.over.chain \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg19_lifted.bed \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg18_liftFail.bed
  cat ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg19_native.bed \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg19_lifted.bed | \
  sed 's/chr//g' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg19_merged.bed
done
#Write to raw CNV files
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/Talkowski_GERM.${CNV}.raw.bed
  sort -Vk1,1 -k2,2n -k3,3n ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.${CNV}.hg19_merged.bed | \
  awk -v OFS="\t" -v CNV=${CNV} '{ print $1, $2, $3, "Talkowski_GERM_"CNV"_"NR, CNV, $6, "22521361" }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/Talkowski_GERM.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/*

#####Assign HPO terms per patient
#Make master list of HPO term matches
cat ${WRKDIR}/data/HPO_map/HPO_map.tsv ${WRKDIR}/bin/rCNVmap/misc/supplementary_HPO_query_terms.txt | \
awk '{ if ($3!=".") print $0 }' > ${WRKDIR}/data/HPO_map/HPO_map.modified.tsv
#Make master list of all patients & phenotypes
cat ${TMPDIR}/../misc_CNVs/TGDB_CNVs.CHB_GeneDX.phenotypes.list \
${TMPDIR}/../misc_CNVs/Coe_Cooper_case_phenotypes.list \
<( printf 'PGC_UNK\tSCHIZOPHRENIA\n%.0s' {1..21094} | awk -v OFS="\t" '{ print $1"_"NR, $2 }' ) \
<( paste <( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | cut -f1 | sort | uniq ) \
<( sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | cut -f1 | sort | uniq | awk -v FS="-" '{ print $2 }' | \
  sed -f ${WRKDIR}/bin/rCNVmap/misc/TCGA_TSS_linkers.sed | awk '{ print toupper($1) }' ) ) > \
${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.list
# <( printf 'Huang_UNK\tTOURETTE\n%.0s' {1..2434} | awk -v OFS="\t" '{ print $1"_"NR, $2 }' ) \
#Iteratively assign HPO terms to each sample
if [ -e ${TMPDIR}/ID_HPO_links.tmp ]; then
  rm ${TMPDIR}/ID_HPO_links.tmp
fi
while read term old HPO count; do
  echo ${term}
  fgrep ${term} ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.list | cut -f1 | \
  awk -v OFS="\t" -v HPO=${HPO} '{ print $1, HPO }' >> ${TMPDIR}/ID_HPO_links.tmp
done < ${WRKDIR}/data/HPO_map/HPO_map.modified.tsv
while read SID PHENOS; do
  HPO=$( fgrep -w ${SID} ${TMPDIR}/ID_HPO_links.tmp | cut -f2 | sed 's/,/\n/g' | sort | uniq | paste -s -d, )
  if [ -z ${HPO} ]; then
    HPO="0000000"
  fi
  echo -e "${SID}\t${HPO}\t${PHENOS}"
done < ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.list > \
${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list
#Get breakdown of number of terms per patient
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | \
sed 's/,/\t/g' | awk '{ print NF }' | sort | uniq -c | sort -nk2,2
#Get breakdown of number of patients per term
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | \
sed 's/,/\n/g' | sort | uniq -c | sort -nrk1,1 | awk -v OFS="\t" '{ print $1, $2 }'
#Get breakdown of number of patients per term (no cancer)
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | \
fgrep -v "0002664" | sed 's/,/\n/g' | sort | uniq -c | sort -nrk1,1 | awk -v OFS="\t" '{ print $1, $2 }'

#####Get counts of patients per analysis group
while read group eti tier description include exclude color n; do
  echo -e "${group}"
  if [ ${eti} == "CNCR" ]; then
    echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
    <( fgrep TCGA ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | cut -f2 ) | \
    fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
  else
    echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
    <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
    fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
  fi
  echo -e "${description}"
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | paste - - -

#####Assign HPO terms to all CNVs
#Cut HPO mappings 
cut -f1,3 ${WRKDIR}/data/HPO_map/HPO_map.modified.tsv > \
${TMPDIR}/HPO_mappings_for_application.txt
#Create list of HPO mappings for all germline CNVs
for cohort in SSC Coe Talkowski PGC Huang; do
  for CNV in DEL DUP; do
    #Cut vector of phenotypes
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | \
    fgrep -v "#" | cut -f6 | tr -d "\.\:\/\\\(\)\'\<\>\~\!\@\#\$\\%\^\&\*\{\}\+\?\,\[\]\-\=" > \
    ${TMPDIR}/${cohort}_${CNV}_GERM.phenos.list
    #Launch HPOmapper
    bsub -q normal -sla miket_sc -J ${cohort}_${CNV}_GERM_HPOmapper \
    "${WRKDIR}/bin/rCNVmap/bin/HPOmapper_helper.R \
    ${TMPDIR}/${cohort}_${CNV}_GERM.phenos.list \
    ${TMPDIR}/HPO_mappings_for_application.txt \
    ${TMPDIR}/${cohort}_${CNV}_GERM.HPOs.list"
  done
done
#Create list of HPO mappings for all cancer CNVs
for cohort in TCGA; do
  for CNV in DEL DUP; do
    #Cut vector of phenotypes
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | \
    fgrep -v "#" | cut -f6 | tr -d "\.\:\/\\\(\)\'\<\>\~\!\@\#\$\\%\^\&\*\{\}\+\?\,\[\]\-\=" > \
    ${TMPDIR}/${cohort}_${CNV}_CNCR.phenos.list
    #Launch HPOmapper
    bsub -q normal -sla miket_sc -J ${cohort}_${CNV}_CNCR_HPOmapper \
    "${WRKDIR}/bin/rCNVmap/bin/HPOmapper_helper.R \
    ${TMPDIR}/${cohort}_${CNV}_CNCR.phenos.list \
    ${TMPDIR}/HPO_mappings_for_application.txt \
    ${TMPDIR}/${cohort}_${CNV}_CNCR.HPOs.list"
  done
done
#Sanity-check HPO mapping output
for cohort in SSC Coe Talkowski PGC Huang; do
  for CNV in DEL DUP; do
    echo -e "${cohort}\t${CNV}"
    #Print number of columns in raw CNV file
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | \
    awk '{ print NF }' | head -n2 | tail -n1
    #Print length of total CNV df
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | \
    fgrep -v "#" | wc -l
    #Print length of phenotypes vector
    cat ${TMPDIR}/${cohort}_${CNV}_GERM.phenos.list | wc -l
    #Print length of HPO list with two columns
    zcat ${TMPDIR}/${cohort}_${CNV}_GERM.HPOs.list.gz | fgrep -v "#" | \
    awk '{ if (NF==2) print $0 }' | wc -l
  done | paste - - - - -
done
for cohort in TCGA; do
  for CNV in DEL DUP; do
    echo -e "${cohort}\t${CNV}"
    #Print number of columns in raw CNV file
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | \
    awk '{ print NF }' | head -n2 | tail -n1
    #Print length of total CNV df
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | \
    fgrep -v "#" | wc -l
    #Print length of phenotypes vector
    cat ${TMPDIR}/${cohort}_${CNV}_CNCR.phenos.list | wc -l
    #Print length of HPO list with two columns
    zcat ${TMPDIR}/${cohort}_${CNV}_CNCR.HPOs.list.gz | fgrep -v "#" | \
    awk '{ if (NF==2) print $0 }' | wc -l
  done | paste - - - - -
done
#Assign HPO mappings to CNV files
for cohort in SSC Coe Talkowski PGC Huang; do
  for CNV in DEL DUP; do
    cat \
    <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | head -n1 ) \
    <( paste <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | \
                fgrep -v "#" | cut -f1-5 ) \
             <( zcat ${TMPDIR}/${cohort}_${CNV}_GERM.HPOs.list.gz | fgrep -v "#" | cut -f1 ) \
             <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | \
                fgrep -v "#" | cut -f7 ) ) | sed 's/,/;/g' > \
    ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed
    gzip -f ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed
  done
done
for cohort in TCGA; do
  for CNV in DEL DUP; do
    cat \
    <( zcat ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_GERM.${CNV}.raw.bed.gz | head -n1 ) \
    <( paste <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | \
                fgrep -v "#" | cut -f1-5 ) \
             <( zcat ${TMPDIR}/${cohort}_${CNV}_CNCR.HPOs.list.gz | fgrep -v "#" | cut -f1 ) \
             <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | \
                fgrep -v "#" | cut -f7 ) ) | sed 's/,/;/g' > \
    ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed
    gzip -f ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed
  done
done

#####Merge all germline CNVs and run bedcluster
# if [ -e ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV ]; then
#   rm -r ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
# fi
# mkdir ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
for CNV in DEL DUP; do
  #Create master list of all CNVs
  # for cohort in Shaikh Suktitipat Uddin Vogler Cooper SSC Coe Talkowski PGC TCGA; do
  #New (3/1/17) Exclude Coe and TCGA control CNVs due to being substantial CNV size + count outliers
  for cohort in Shaikh Suktitipat Uddin Vogler Cooper SSC Talkowski PGC Coe; do
    if [ -e ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz ]; then
      zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | fgrep -v "#" | \
      awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
    fi
    if [ ${cohort} != "Coe" ]; then
      if [ -e ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CTRL.${CNV}.raw.bed.gz ]; then
        zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CTRL.${CNV}.raw.bed.gz | fgrep -v "#" | \
        awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
      fi
    fi
  done | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.bed
  #Run bedtools intersect (50% recip) and require both breakpoints ±20kb
  bedtools intersect -r -f 0.5 -wa -wb \
  -a ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.bed \
  -b ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.bed | \
  awk -v d=20000 '{ if ($2-$8<=d && $2-$8>=-d && $3-$9<=d && $3-$9>=-d) print $0 }' > \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.all_vs_all.bed
  #Run bedcluster
  cut -f4 ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.VIDs.list
  bsub -q big -R 'rusage[mem=50000]' -M 50000 -sla miket_sc -J all_germline_${CNV}_merge \
  "source /apps/lab/miket/anaconda/4.0.5/envs/collins_py3/bin/activate collins_py3; \
  /data/talkowski/rlc47/code/svcf/scripts/bedcluster -p all_germline_${CNV} -m \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.VIDs.list \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.all_vs_all.bed \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed"
done

#####Run bedcluster on cancer CNVs
for CNV in DEL DUP; do
  #Create master list of all CNVs
  for cohort in TCGA; do
    if [ -e ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz ]; then
      zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CNCR.${CNV}.raw.bed.gz | fgrep -v "#" | \
      awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
    fi
  done | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed
  #Run bedtools intersect (50% recip) and require both breakpoints ±20kb
  bedtools intersect -r -f 0.5 -wa -wb \
  -a ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed \
  -b ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed | \
  awk -v d=20000 '{ if ($2-$8<=d && $2-$8>=-d && $3-$9<=d && $3-$9>=-d) print $0 }' > \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.all_vs_all.bed
  #Run bedcluster
  cut -f4 ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed | \
  sort | uniq > ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.VIDs.list
  bsub -q big -R 'rusage[mem=50000]' -M 50000 -sla miket_sc -J all_cancer_${CNV}_merge \
  "source /apps/lab/miket/anaconda/4.0.5/envs/collins_py3/bin/activate collins_py3; \
  /data/talkowski/rlc47/code/svcf/scripts/bedcluster -p all_cancer_${CNV} -m \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.VIDs.list \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.all_vs_all.bed \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed"
done

# #####New (3/1/17) Exclude Coe and TCGA control CNVs due to being substantial CNV size + count outliers
# for CNV in DEL DUP; do
#   fgrep -ve "Coe_CTRL" -e "TCGA_CTRL" ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
#   ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.noCoeTCGACTRL.bed
# done



######################################################
# GROUP 1 (E2): VF < 1%, size > 50kb, size < 5Mb
######################################################

#####Filter merged CNVs across cohort (max VF)
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2 ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/
max_VF=0.01
for CNV in DEL DUP; do
  #Germline
  awk -v max_VF=${max_VF} '{ if (($9/104909)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/germline_${CNV}.merged.maxVF.bed
  #Cancer
  cat ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/cancer_${CNV}.merged.maxVF.bed
done

#####Filter merged CNVs per study (VF < 10% per study) & reassign original coordinates
max_VF=0.1
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
    study_base=$( echo "${study}" | cut -f1 -d_ )
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*.maxVF.bed | fgrep "${study}_${CNV}" | \
    cut -f7 | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
    awk -v max_VF=${max_VF} -v n=${n} '{ if (($2/n)<max_VF) print $1 }' | \
    fgrep -wf - ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*.maxVF.bed | \
    cut -f2- -d\: | cut -f4 | fgrep "${study}_${CNV}" | sort | uniq | \
    fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  done < <( fgrep -ve "TCGA" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/cancer_${CNV}.merged.maxVF.bed | \
  cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed.gz ) >> \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
done

#####Filter merged CNVs on minimum size
min_size=50000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
    awk -v min_size=${min_size} '{ if ($3-$2>min_size) print $0 }' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed | \
    bedtools intersect -v -f 0.5 -a - -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on exclusion loci
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
    sed -e 's/^x/X/g' -e 's/^y/Y/g' -e 's/^MT/M/g' -e 's/^5_/5/g' -e 's/^16_/16/g' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed | \
    bedtools intersect -v -f 0.5 -a - \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
    fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
    awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Make master germline, control, and cancer files (post-filter)
for pheno in GERM CTRL CNCR; do
  #Require CNCR samples to come from TCGA only
  if [ ${pheno} == "CNCR" ]; then
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  else
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  fi    
done
#Gzip all filtered CNV files
gzip ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/*

#####Final merger per analysis group, filtered on max size & unfiltered on max size
if [ -e ${WRKDIR}/data/CNV/CNV_MASTER/ ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_MASTER/
fi
mkdir ${WRKDIR}/data/CNV/CNV_MASTER/
while read group eti tier descrip include exclude color n; do
  echo ${group}
  if [ -e ${WRKDIR}/data/CNV/CNV_MASTER/${group} ]; then
    rm -rf ${WRKDIR}/data/CNV/CNV_MASTER/${group}
  fi
  mkdir ${WRKDIR}/data/CNV/CNV_MASTER/${group}
  if [ ${exclude} != "NA" ]; then
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E2.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E2.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E2.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E2.GRCh37.all.bed
    done
  else
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E2.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E2.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E2.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E2.GRCh37.all.bed
    done
  fi
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Make merged CNV (DEL+DUP) set
while read group eti tier descrip include exclude color n; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.E2.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.E2.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E2.GRCh37.all.bed
  #No max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E2.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.noMaxSize.E2.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.noMaxSize.E2.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E2.GRCh37.all.bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Generate coding, noncoding, intergenic, and haplosufficient CNV callsets
while read group eti tier descrip include exclude color n; do
  echo ${group}
  for class in CNV DEL DUP; do
    echo ${class}
    #Coding, with max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.coding.bed
    #Noncoding, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.noncoding.bed
    #Intergenic, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.intergenic.bed
    #Haplosufficient, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E2.GRCh37.haplosufficient.bed
    #Coding, no max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.coding.bed
    #Noncoding, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.noncoding.bed
    #Intergenic, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.intergenic.bed
    #Haplosufficient, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E2.GRCh37.haplosufficient.bed
  done
  gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/${group}/*E2*bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )



######################################################
# GROUP 2 (E3): VF < 0.1%, size > 50kb, size < 5Mb
######################################################

#####Filter merged CNVs across cohort (max VF)
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3 ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/
max_VF=0.001
for CNV in DEL DUP; do
  #Germline
  awk -v max_VF=${max_VF} '{ if (($9/104909)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/germline_${CNV}.merged.maxVF.bed
  #Cancer
  cat ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/cancer_${CNV}.merged.maxVF.bed
done

#####Filter merged CNVs per study (VF < 1% per study) & reassign original coordinates
max_VF=0.01
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
    study_base=$( echo "${study}" | cut -f1 -d_ )
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*.maxVF.bed | fgrep "${study}_${CNV}" | \
    cut -f7 | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
    awk -v max_VF=${max_VF} -v n=${n} '{ if (($2/n)<max_VF) print $1 }' | \
    fgrep -wf - ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*.maxVF.bed | \
    cut -f2- -d\: | cut -f4 | fgrep "${study}_${CNV}" | sort | uniq | \
    fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  done < <( fgrep -ve "TCGA" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/cancer_${CNV}.merged.maxVF.bed | \
  cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed.gz ) >> \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
done

#####Filter merged CNVs on minimum size
min_size=50000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
    awk -v min_size=${min_size} '{ if ($3-$2>min_size) print $0 }' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed | \
    bedtools intersect -v -f 0.5 -a - -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on exclusion loci
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
    sed -e 's/^x/X/g' -e 's/^y/Y/g' -e 's/^MT/M/g' -e 's/^5_/5/g' -e 's/^16_/16/g' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed | \
    bedtools intersect -v -f 0.5 -a - \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
    fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
    awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Make master germline, control, and cancer files (post-filter)
for pheno in GERM CTRL CNCR; do
  #Require CNCR samples to come from TCGA only
  if [ ${pheno} == "CNCR" ]; then
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  else
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  fi    
done
#Gzip all filtered CNV files
gzip ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/*

#####Final merger per analysis group, filtered on max size & unfiltered on max size
while read group eti tier descrip include exclude color n; do
  echo ${group}
  if [ ${exclude} != "NA" ]; then
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E3.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E3.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E3.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E3.GRCh37.all.bed
    done
  else
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E3.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E3.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E3.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E3/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E3.GRCh37.all.bed
    done
  fi
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Make merged CNV (DEL+DUP) set
while read group eti tier descrip include exclude color n; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E3.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.E3.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.E3.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E3.GRCh37.all.bed
  #No max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E3.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.noMaxSize.E3.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.noMaxSize.E3.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E3.GRCh37.all.bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Generate coding, noncoding, intergenic, and haplosufficient CNV callsets
while read group eti tier descrip include exclude color n; do
  echo ${group}
  for class in CNV DEL DUP; do
    echo ${class}
    #Coding, with max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.coding.bed
    #Noncoding, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.noncoding.bed
    #Intergenic, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.intergenic.bed
    #Haplosufficient, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E3.GRCh37.haplosufficient.bed
    #Coding, no max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.coding.bed
    #Noncoding, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.noncoding.bed
    #Intergenic, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.intergenic.bed
    #Haplosufficient, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E3.GRCh37.haplosufficient.bed
  done
  gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/${group}/*E3*bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )









######################################################
# GROUP 3 (E4): VF < 0.01%, size > 50kb, size < 5Mb
######################################################

#####Filter merged CNVs across cohort (max VF)
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4 ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/
max_VF=0.0001
for CNV in DEL DUP; do
  #Germline
  awk -v max_VF=${max_VF} '{ if (($9/104909)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/germline_${CNV}.merged.maxVF.bed
  #Cancer
  cat ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/cancer_${CNV}.merged.maxVF.bed
done

#####Filter merged CNVs per study (VF < 0.01% or n=1 per study) & reassign original coordinates
max_VF=0.001
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
    study_base=$( echo "${study}" | cut -f1 -d_ )
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*.maxVF.bed | fgrep "${study}_${CNV}" | \
    cut -f7 | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
    awk -v max_VF=${max_VF} -v n=${n} '{ if (($2/n)<max_VF || $2==1) print $1 }' | \
    fgrep -wf - ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*.maxVF.bed | \
    cut -f2- -d\: | cut -f4 | fgrep "${study}_${CNV}" | sort | uniq | \
    fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  done < <( fgrep -ve "TCGA" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/cancer_${CNV}.merged.maxVF.bed | \
  cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed.gz ) >> \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
done

#####Filter merged CNVs on minimum size
min_size=50000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
    awk -v min_size=${min_size} '{ if ($3-$2>min_size) print $0 }' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed | \
    bedtools intersect -v -f 0.5 -a - -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on exclusion loci
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
    sed -e 's/^x/X/g' -e 's/^y/Y/g' -e 's/^MT/M/g' -e 's/^5_/5/g' -e 's/^16_/16/g' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed | \
    bedtools intersect -v -f 0.5 -a - \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
    fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
    awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
done

#####Make master germline, control, and cancer files (post-filter)
for pheno in GERM CTRL CNCR; do
  #Require CNCR samples to come from TCGA only
  if [ ${pheno} == "CNCR" ]; then
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  else
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters.bed
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    done
  fi    
done
#Gzip all filtered CNV files
gzip ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/*

#####Final merger per analysis group, filtered on max size & unfiltered on max size
while read group eti tier descrip include exclude color n; do
  echo ${group}
  if [ ${exclude} != "NA" ]; then
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E4.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E4.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E4.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E4.GRCh37.all.bed
    done
  else
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E4.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E4.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E4.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E4/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E4.GRCh37.all.bed
    done
  fi
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Make merged CNV (DEL+DUP) set
while read group eti tier descrip include exclude color n; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E4.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.E4.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.E4.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.E4.GRCh37.all.bed
  #No max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E4.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.noMaxSize.E4.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.noMaxSize.E4.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E4.GRCh37.all.bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Generate coding, noncoding, intergenic, and haplosufficient CNV callsets
while read group eti tier descrip include exclude color n; do
  echo ${group}
  for class in CNV DEL DUP; do
    echo ${class}
    #Coding, with max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.coding.bed
    #Noncoding, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.noncoding.bed
    #Intergenic, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.intergenic.bed
    #Haplosufficient, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.E4.GRCh37.haplosufficient.bed
    #Coding, no max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.coding.bed
    #Noncoding, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.noncoding.bed
    #Intergenic, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.intergenic.bed
    #Haplosufficient, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.E4.GRCh37.haplosufficient.bed
  done
  gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/${group}/*E4*bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )









######################################################
# GROUP 4 (N1): n=1, size > 50kb, size < 5Mb
######################################################

if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/
#####Identify N1s
#Apply 50kb min size cutoff, apply blacklist filter
min_size=50000
for CNV in DEL DUP; do
  awk -v min_size=${min_size} '{ if ($9==1 && $3-$2>min_size) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
  bedtools intersect -v -f 0.5 -a - \
  -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.noMaxSize.bed
done
#Supplement with SSC N1s that appear in up to 2 samples from the SSC (and no other studies)
for CNV in DEL DUP; do
  while read VID; do
    others=$( fgrep -w ${VID} ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
              fgrep -v SSC | wc -l )
    if [ ${others} == 0 ]; then
      fgrep -w ${VID} ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed
    fi
  done < <( fgrep SSC ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
    bedtools intersect -v -f 0.5 -a - -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | \
  awk -v min_size=${min_size} '{ if ($9==2 && $3-$2>min_size) print $7 }' | sort | uniq ) >> \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.noMaxSize.bed
done

#Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.noMaxSize.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.bed
done
#Convert to per-study files
for CNV in DEL DUP; do
  while read study n PMID; do
    study_base=$( echo "${study}" | cut -f1 -d_ )
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    #No max size
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/${study}_${CNV}.noMaxSize.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.noMaxSize.bed | fgrep "${study}_${CNV}" | \
    cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/${study}_${CNV}.noMaxSize.bed
    #With max size
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/${study}_${CNV}.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/germline_${CNV}.bed | fgrep "${study}_${CNV}" | \
    cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/${study}_${CNV}.bed
  done < <( fgrep -ve "TCGA_CTRL" -e "Coe_CTRL" ${WRKDIR}/lists/Studies_SampleSizes.list )
  #CNCR, no max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/TCGA_CNCR_${CNV}.noMaxSize.bed
  zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_CNCR_${CNV}.final_filters_noMaxSize.bed.gz | \
  fgrep -v "#" >> ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/TCGA_CNCR_${CNV}.noMaxSize.bed
  #CNCR, with max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/TCGA_CNCR_${CNV}.bed
  zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_E2/ALL_CNCR_${CNV}.final_filters.bed.gz | \
  fgrep -v "#" >> ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/TCGA_CNCR_${CNV}.bed
done
#Make master N1 germline, control, and cancer files (post-filter)
for pheno in GERM CTRL CNCR; do
  #Require CNCR samples to come from TCGA only
  if [ ${pheno} == "CNCR" ]; then
    for CNV in DEL DUP; do
      #No max size
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/*${pheno}_${CNV}.noMaxSize.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      #With max size
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/*${pheno}_${CNV}.bed | \
      fgrep -v "#" | fgrep TCGA | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters.bed
    done
  else
    for CNV in DEL DUP; do
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/*${pheno}_${CNV}.noMaxSize.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
      #With max size
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters.bed
      cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/*${pheno}_${CNV}.bed | \
      fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${pheno}_${CNV}.final_filters.bed
    done
  fi    
done
#Gzip all
gzip ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/*

#####Final N1 merger per analysis group, filtered on max size & unfiltered on max size
while read group eti tier descrip include exclude color n; do
  echo ${group}
  if [ ${exclude} != "NA" ]; then
    for CNV in DEL DUP; do
      #Filtered on max size (main N1 set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.N1.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.N1.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.N1.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.N1.GRCh37.all.bed
    done
  else
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set for analysis)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.N1.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.N1.GRCh37.all.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.N1.GRCh37.all.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_N1/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.N1.GRCh37.all.bed
    done
  fi
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Make merged CNV (DEL+DUP) set
while read group eti tier descrip include exclude color n; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.N1.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.N1.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.N1.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.N1.GRCh37.all.bed
  #No max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.N1.GRCh37.all.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DEL.noMaxSize.N1.GRCh37.all.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.DUP.noMaxSize.N1.GRCh37.all.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.N1.GRCh37.all.bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Generate coding, noncoding, intergenic, and haplosufficient CNV callsets
while read group eti tier descrip include exclude color n; do
  echo ${group}
  for class in CNV DEL DUP; do
    echo ${class}
    #Coding, with max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.coding.bed
    #Noncoding, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.noncoding.bed
    #Intergenic, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.intergenic.bed
    #Haplosufficient, with max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.N1.GRCh37.haplosufficient.bed
    #Coding, no max size
    bedtools intersect -u -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.coding.bed
    #Noncoding, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.noncoding.bed
    #Intergenic, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.intergenic.bed
    #Haplosufficient, no max size
    bedtools intersect -v -wa \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.notHaplosufficient.bed | \
    sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.all.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${class}.noMaxSize.N1.GRCh37.haplosufficient.bed
  done
  gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/${group}/*N1*bed
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )