#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Ryan's working code for analyses & processing conducted on ERISone (PHC HPC cluster)

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

#####Create directory tree
mkdir ${WRKDIR}/lists
mkdir ${WRKDIR}/bin
mkdir ${WRKDIR}/data
mkdir ${WRKDIR}/data/CNV
mkdir ${WRKDIR}/data/CNV/CNV_RAW
mkdir ${WRKDIR}/data/CNV/CNV_MASTER
mkdir ${WRKDIR}/data/CNV/CNV_DENSITY
mkdir ${WRKDIR}/data/plot_data
mkdir ${WRKDIR}/data/annotations
mkdir ${WRKDIR}/data/annotations/exonExclusion
mkdir ${WRKDIR}/analysis
mkdir ${WRKDIR}/analysis/BIN_CNV_pileups
mkdir ${WRKDIR}/analysis/BIN_CNV_burdens
mkdir ${WRKDIR}/analysis/BIN_CNV_permutation
mkdir ${WRKDIR}/analysis/Final_Loci
mkdir ${WRKDIR}/analysis/Final_Loci/significant
mkdir ${WRKDIR}/analysis/Functional_Enrichments
mkdir ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion

#####Gather raw SSC CNVs
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
#Split by case/control
for CNV in DEL DUP; do
  #Controls
  cat <( echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" ) \
  <( awk '$4 ~ /fa|mo|s/ { print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | sed 's/\?//g' | \
  awk -v CNV=${CNV} -v OFS="\t" '{ print $1, $2, $3, "SSC_CTRL_"CNV"_"NR, CNV, "CTRL", "26402605" }' ) > \
  ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_CTRL.${CNV}.raw.bed
  #Cases
  cat <( echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" ) \
  <( awk '$4 ~ /p/ { print $0 }' ${TMPDIR}/SSC_CNVs.p10E9.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | sed 's/\?//g' | \
  awk -v CNV=${CNV} -v OFS="\t" '{ print $1, $2, $3, "SSC_DD_"CNV"_"NR, CNV, "DD", "26402605" }' ) > \
  ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/SSC_DD.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/SSC_CNVs/*.bed

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

#####Gather Itsara CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/Itsara_CTRL.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tCTRL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Itsara_CNVs/Itsara.DUP.raw.bed.gz ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Itsara_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/Itsara_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Itsara_CNVs/*.raw.bed

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
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tCTRL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Shaikh_CNVs/Shaikh.DUP.raw.bed.gz ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Shaikh_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs/Shaikh_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Shaikh_CNVs/*.raw.bed

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
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tCTRL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Suktitipat_Thai_CNVs/Suktitipat.DUP.raw.bed.gz ) | \
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
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tCTRL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Uddin_Ontario_CNVs/Uddin.DUP.raw.bed.gz ) | \
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
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tCTRL\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Vogler_CNVs/Vogler.DUP.raw.bed.gz ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Vogler_CTRL_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs/Vogler_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Vogler_CNVs/*.raw.bed

#####Gather Coe CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe_DD.${CNV}.raw.bed
  while read chr start end sID skip n total PMID; do
    for times in $( seq 1 ${n} ); do
      echo -e "${chr}\t${start}\t${end}\tID\t${CNV}\tDD\t${PMID}"
    done
  done < <( zcat /scratch/miket/rlc47temp/DGV_CNVs_all/Coe_et_al_2014_CNVs/Coe.DUP.raw.bed.gz ) | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "Coe_DD_"CNV"_"NR, $5, $6 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/Coe_DD.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Coe_CNVs/*.raw.bed

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
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_SCZ.${CNV}.raw.bed
  fgrep -w Case ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "PGC_SCZ_"CNV"_"NR, CNV, "SCZ", "27869829" }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_SCZ.${CNV}.raw.bed
  #Control
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_CTRL.${CNV}.raw.bed
  fgrep -w Control ${TMPDIR}/../misc_CNVs/PGC_SCZ_41K_CNV_EGA.sorted.hg19.${CNV}.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v CNV=${CNV} \
  '{ print $1, $2, $3, "PGC_CTRL_"CNV"_"NR, CNV, "CTRL", "27869829" }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/PGC_CTRL.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/PGC_CNVs/*.raw.bed

#####Gather TCGA CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs
#Filter CNV classes by log2 ratios
for CNV in DEL DUP; do
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed
done
sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF<=-1) print $0 }' | \
sort -Vk2,2 -k3,3 -k4,4n | awk -v OFS="\t" '{ print $2, $3, $4, "TCGA_CNCR_DEL_"NR, "DEL", "CNCR", "24071852" }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.DEL.raw.bed
sed '1d' ${TMPDIR}/../misc_CNVs/all_cancers.seg | awk '{ if ($NF>=0.5849625) print $0 }' | \
sort -Vk2,2 -k3,3 -k4,4n | awk -v OFS="\t" '{ print $2, $3, $4, "TCGA_CNCR_DUP_"NR, "DUP", "CNCR", "24071852" }' >> \
${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.DUP.raw.bed
gzip ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/*.raw.bed

#####Gather Talkowski CNVs
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs
#Download from database (MySQL code below)
# SELECT `Chr`, `Start`, `Stop`, `Diag_state`, `ID`, `Source`, `case_type`, `Ref` \
# FROM all_cnvs WHERE Source = 'CHB' OR Source = 'Cooper_etal' OR Source = 'GENEDX';
#Split by CNV class
sed 's/\ /_/g' ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed2
mv ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed2 \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed
awk -v OFS="\t" -v FS="\t" '{ if ($4=="0,200,0" || $4=="COPY_LOSS" || \
  $4=="COPY_LOSS." || $4=="COPY_LOSSETION" || $4=="COPY_LOSS_(MOS)" || \
  $4=="kiss" || $4=="los" || $4=="loss" || $4=="loss_" || $4=="Loss" || \
  $4=="LOSS" || $4=="MOS_COPY_LOSS" || $4=="NULL_COPY_LOSS") print $0 }' \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.DEL.bed
awk -v OFS="\t" -v FS="\t" '{ if ($4=="220,0,0" || $4=="255,0,0" || $4=="AMP" || \
  $4=="AMPLIFICATION" || $4=="COPY_GAIN" || $4=="COPY_GAIN_" || \
  $4=="COPY_GAIN_(MOS?)" || $4=="COPY_GAIN_(MOS)" || $4=="dup" || $4=="g" || \
  $4=="gain" || $4=="gain_" || $4=="Gain" || $4=="GAIN" || \
  $4=="MOSAIC_COPY_GAIN" || $4=="MOS_COPY_GAIN") print $0 }' \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.bed > \
${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.DUP.bed
#Split by native reference
for CNV in DEL DUP; do
  for ref in hg18 hg19; do
    awk -v ref=${ref} '{ if ($NF==ref) print $0 }' \
    ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.bed > \
    ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.${ref}_native.bed
  done
done
#Liftover to hg19
for CNV in DEL DUP; do
  liftOver -minMatch=0.5 -bedPlus=4 \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg18_native.bed \
  /data/talkowski/rlc47/src/hg18ToHg19.over.chain \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_lifted.bed \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg18_liftFail.bed
  cat ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_native.bed \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_lifted.bed | \
  sed 's/chr//g' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_merged.bed
done
#Split by case/control
for CNV in DEL DUP; do
  for pheno in CTRL DD; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/Talkowski_${pheno}.${CNV}.raw.bed
  done
  #Controls
  fgrep Cooper_etal ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_merged.bed | \
  fgrep -v WTCCC | awk -v OFS="\t" -v CNV=${CNV} '{ print $1, $2, $3, "Talkowski_CTRL_"CNV"_"NR, CNV, "CTRL", 22521361 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/Talkowski_CTRL.${CNV}.raw.bed
  #Cases
  fgrep -v Cooper_etal ${TMPDIR}/../misc_CNVs/TGDB_CNVs.noExcluded.CHB_Cooper_GeneDX.${CNV}.hg19_merged.bed | \
  awk -v OFS="\t" -v CNV=${CNV} '{ print $1, $2, $3, "Talkowski_DD_"CNV"_"NR, CNV, "DD", 22521361 }' >> \
  ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/Talkowski_DD.${CNV}.raw.bed
done
gzip ${WRKDIR}/data/CNV/CNV_RAW/Talkowski_CNVs/*

#####Merge all germline CNVs and run bedcluster
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
for CNV in DEL DUP; do
  #Create master list of all CNVs
  for cohort in SSC Itsara Shaikh Suktitipat Uddin Vogler Coe PGC Talkowski; do
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/*.${CNV}.raw.bed.gz | fgrep -v "#" | \
    awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
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
  bsub -q normal -sla miket_sc -J all_germline_${CNV}_merge \
  "/data/talkowski/rlc47/code/svcf/scripts/bedcluster -p all_germline_${CNV} -m \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.VIDs.list \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.pre_merge.all_vs_all.bed \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed"
done

#####Run bedcluster on cancer CNVs
for CNV in DEL DUP; do
  #Create master list of all CNVs
  for cohort in TCGA; do
    zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/*.${CNV}.raw.bed.gz | fgrep -v "#" | \
    awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
  done | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed
  #Run bedtools intersect (50% recip) and require both breakpoints ±20kb
  bedtools intersect -r -f 0.5 -wa -wb \
  -a ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed \
  -b ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed | \
  awk -v d=20000 '{ if ($2-$8<=d && $2-$8>=-d && $3-$9<=d && $3-$9>=-d) print $0 }' > \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.all_vs_all.bed
  #Run bedcluster
  cut -f4 ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.VIDs.list
  bsub -q normal -sla miket_sc -J all_cancer_${CNV}_merge \
  "/data/talkowski/rlc47/code/svcf/scripts/bedcluster -p all_cancer_${CNV} -m \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.VIDs.list \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.pre_merge.all_vs_all.bed \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed"
done

#####Filter merged CNVs across cohort (min size & max VF)
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/
min_size=20000
max_VF=0.01
for CNV in DEL DUP; do
  #Germline
  awk -v min_size=${min_size} -v max_VF=${max_VF} \
  '{ if ($3-$2>min_size && ($9/112053)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/germline_${CNV}.merged.minSize_maxVF.bed
  #Cancer
  awk -v min_size=${min_size} -v max_VF=${max_VF} \
  '{ if ($3-$2>min_size && ($9/10844)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/cancer_${CNV}.merged.minSize_maxVF.bed
done

#####Filter merged CNVs per study (max VF per study)
max_VF=0.01
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.bed
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*.minSize_maxVF.bed | fgrep "${study}_${CNV}" | \
    cut -f7 | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
    awk -v max_VF=${max_VF} -v n=${n} '{ if (($2/n)<max_VF) print $1 }' | \
    fgrep -wf - <( cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*.minSize_maxVF.bed ) | \
    fgrep "${study}_${CNV}" | awk -v OFS="\t" -v PMID=${PMID} -v CNV=${CNV} -v pheno=${pheno} \
    '{ print $1, $2, $3, $4, CNV, pheno, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Filter merged CNVs on exclusion loci
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.bed
    bedtools intersect -v -f 0.5 \
    -a ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.bed \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.maxSize.bed
    fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.bed | \
    awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.maxSize.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Final merger per cohort, filtered on max size & unfiltered on max size
if [ -e ${WRKDIR}/data/CNV/CNV_MASTER/ ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_MASTER/
fi
mkdir ${WRKDIR}/data/CNV/CNV_MASTER/
for group in DD SCZ CNCR CTRL; do
  for CNV in DEL DUP; do
    #Filtered on max size (main CNV set)
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*_${group}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.maxSize.bed | \
    fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' >> \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
    #Unfiltered on max size (used for size distribs)
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*_${group}_${CNV}.merged.minSize_maxVF.maxVF.blacklisted.bed | \
    fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' >> \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
  done
done

#####Make merged CNV (DEL+DUP) set
for group in DD SCZ CNCR CTRL; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.CNV.GRCh37.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.DEL.GRCh37.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.DUP.GRCh37.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.CNV.GRCh37.bed
  #No max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.CNV.noMaxSize.GRCh37.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.DEL.noMaxSize.GRCh37.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.DUP.noMaxSize.GRCh37.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/${group}.CNV.noMaxSize.GRCh37.bed
done

#####Make merged DD+SCZ set
for CNV in CNV DEL DUP; do
  #With max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/DD_SCZ.${CNV}.GRCh37.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/DD.${CNV}.GRCh37.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/SCZ.${CNV}.GRCh37.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/DD_SCZ.${CNV}.GRCh37.bed
  #Without max size
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_MASTER/DD_SCZ.${CNV}.noMaxSize.GRCh37.bed
  cat ${WRKDIR}/data/CNV/CNV_MASTER/DD.${CNV}.noMaxSize.GRCh37.bed \
  ${WRKDIR}/data/CNV/CNV_MASTER/SCZ.${CNV}.noMaxSize.GRCh37.bed | fgrep -v "#" | \
  sort -Vk1,1 -k2,2n -k3,3n >> \
  ${WRKDIR}/data/CNV/CNV_MASTER/DD_SCZ.${CNV}.noMaxSize.GRCh37.bed
done

#####Generate coding and noncoding CNV callsets
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for class in CNV DEL DUP; do
    #Coding, with max size
    awk '$4 !~ /\-AS1/ { print $0 }' ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed | \
    bedtools intersect -u -wa -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.bed | sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.coding.bed
    #Noncoding, with max size
    awk '$4 !~ /\-AS1/ { print $0 }' ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed | \
    bedtools intersect -v -wa -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.bed | sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.GRCh37.noncoding.bed
    #Coding, no max size
    awk '$4 !~ /\-AS1/ { print $0 }' ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed | \
    bedtools intersect -u -wa -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.bed | sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.coding.bed
    #Noncoding, no max size
    awk '$4 !~ /\-AS1/ { print $0 }' ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed | \
    bedtools intersect -v -wa -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.bed | sort -Vk1,1 -k2,2n -k3,3n | \
    cat <( head -n1 ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.bed ) - > \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${class}.noMaxSize.GRCh37.noncoding.bed
  done
done
gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/*bed

#####Generate CNV density per callset
if [ -e ${WRKDIR}/data/CNV/CNV_DENSITY ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_DENSITY
fi
mkdir ${WRKDIR}/data/CNV/CNV_DENSITY
for group in CTRL DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for CNV in DEL DUP CNV; do
    echo ${CNV}
    bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
    -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz > \
    ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
    gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
    bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
    -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz > \
    ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
    gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
    bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
    -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz > \
    ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
    gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
  done
done

#####Sanity-check filtering
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    echo -e "${group}.${CNV}"
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | wc -l
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l
  done | paste - - - - | awk -v OFS="\t" '{ print $0, $3+$4 }'
done

#####Gather distributions of CNV per cohort per phenotype
for CNV in CNV DEL DUP; do
  echo -e "\n${CNV}\n"
  for pheno in CTRL DD SCZ CNCR; do
    while read cohort; do
      echo -e "${cohort}\t${pheno}\t${CNV}"
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    done < <( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | \
      fgrep -v "#" | cut -f4 | cut -f1 -d_ | sort | uniq ) | paste - - - - - - -
    for dummy in 1; do
      echo -e "ALL\t${pheno}\t${CNV}"
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    done | paste - - - - - - -
  done
done
#Get median sizes of all CNVs by class & germline/all
for CNV in CNV DEL DUP; do
  echo -e "\n${CNV}\n"
  for dummy in 1; do
    for pheno in DD SCZ CTRL; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    for pheno in DD SCZ CTRL; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    for pheno in DD SCZ CTRL; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  done | paste - - -
  for dummy in 1; do
    for pheno in DD SCZ CTRL CNCR; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    for pheno in DD SCZ CTRL CNCR; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    for pheno in DD SCZ CTRL CNCR; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
    done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
  done | paste - - -
done

#####Get reverse CDF of CNV sizes by phenotype
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP; do
    for dummy in 1; do
      echo -e "${group}.${CNV}"
      total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | wc -l )
      for suffix in 000 0000 00000; do
        for prefix in $( seq 10 1 90 ); do
          filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | \
          awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
          echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
        done
      done | paste -s
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    done | paste -s
  done
done > ${WRKDIR}/data/plot_data/CNV_size_bySize.all.txt
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP; do
    for dummy in 1; do
      echo -e "${group}.${CNV}"
      total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l )
      for suffix in 000 0000 00000; do
        for prefix in $( seq 10 1 90 ); do
          filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | \
          awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
          echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
        done
      done | paste -s
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    done | paste -s
  done
done > ${WRKDIR}/data/plot_data/CNV_size_bySize.coding.txt
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP; do
    for dummy in 1; do
      echo -e "${group}.${CNV}"
      total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l )
      for suffix in 000 0000 00000; do
        for prefix in $( seq 10 1 90 ); do
          filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
          awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
          echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
        done
      done | paste -s
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
      awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
    done | paste -s
  done
done > ${WRKDIR}/data/plot_data/CNV_size_bySize.noncoding.txt

#####Run TBRden pileup for all tissue types and all phenotypes (coding & noncoding CNVs)
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_coding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.coding.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
  done
done

#####Run TBRden analysis
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  case ${group} in
    DD)
      color="green"
      ;;
    SCZ)
      color="purple"
      ;;
    DD_SCZ)
      color="blue"
      ;;
    CNCR)
      color="orange"
      ;;
  esac
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    #Parallelize analyses (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_all ${color}"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_coding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.coding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.coding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_coding ${color}"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_noncoding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_noncoding ${color}"
  done
done

#####Run 100k CNV shift permutation tests for all comparisons
#Note: initial p-value cutoff used: 0.05/26802 = 1.865532e-06
#This corresponds to the number of non-overlapping autosomal 100kb bins we tested (after blacklisting N-mask, etc)
for group in DD SCZ DD_SCZ CNCR; do
  echo ${group}
  if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split
  for CNV in DEL DUP CNV; do
    echo ${CNV}
    #Coding + noncoding CNVs
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
    #Coding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed
    #Noncoding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
    #Split into partitions of 1k permutations each (x100 per comparison)
    for i in $( seq -w 001 100 ); do
      echo ${i}
      #Coding + noncoding CNVs
      if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i} ]; then
        mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}
      fi
      bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 2 -b 1000000 -N 10 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_all \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
      #Coding CNVs only
      bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 2 -b 100000 -N 45 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_coding \
      -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
      #Noncoding CNVs only
      bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 2 -b 100000 -N 1000 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_noncoding \
      -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
    done
  done
done

#####Collect results from permutation tests
#Collect results across all permutations per group
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      for i in $( seq -w 001 100 ); do
        zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_${filt}.permuted.${i}.bed.gz | \
        fgrep -v "#" | cut -f7 | paste -s
      done > ${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt
      Rscript -e "options(scipen=100);\
      d <- apply(read.table(\"${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt\",header=F),2,sum);\
      write.table(data.frame(\"perms_less_sig\"=100000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/100001),\
        \"${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt\",row.names=F,col.names=T,sep=\"\t\",quote=F)"
      paste <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/001/${group}_vs_CTRL_${CNV}_${filt}.permuted.001.bed.gz | cut -f1-5 ) \
      ${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt > \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed
    done
  done
done
#Compose table of all bins for each comparison
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      cat <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -wb -f 1 -r -a - \
      -b ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed | \
      cut -f1-17,23-25 ) \
      <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -v -f 1 -r -a - \
      -b ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed | \
      awk -v OFS="\t" '{ print $0, "NA\tNA\tNA" }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
      cat <( paste <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
        head -n1 | awk '{ print "#"$0 }' ) <( echo -e "perms_less_sig\tperms_as_or_more_sig\tperm_p" ) ) - > \
      ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
    done
  done
done
#Make master table of all bins with summary stats for each bin & each comparison
for dummy in 1; do
  for group in CTRL DD SCZ CNCR; do
    for CNV in DEL DUP; do
      for filt in coding noncoding; do
        zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
        fgrep -v "#" | cut -f5 > ${TMPDIR}/${group}.${CNV}.${filt}.CNVcounts.tmp
        echo ${TMPDIR}/${group}.${CNV}.${filt}.CNVcounts.tmp
      done
    done
  done
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all coding noncoding; do
        zcat ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed.gz | \
        fgrep -v "#" | cut -f17,20 > ${TMPDIR}/${group}.${CNV}.${filt}.pvals.tmp
        echo ${TMPDIR}/${group}.${CNV}.${filt}.pvals.tmp
      done
    done
  done
done > ${TMPDIR}/list_of_paths.txt
for dummy in 1; do
  zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
  cut -f1-4 | head -n1
  for group in CTRL DD SCZ CNCR; do
    for CNV in DEL DUP; do
      for filt in coding noncoding; do
        echo "${group}.${CNV}.${filt}.count"
      done
    done
  done
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all coding noncoding; do
        echo -e "${group}.${CNV}.${filt}.obs_p\n${group}.${CNV}.${filt}.perm_p"
      done
    done
  done
done | paste -s > ${TMPDIR}/header.txt
paste <( zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
  fgrep -v "#" | cut -f1-4 ) $( cat ${TMPDIR}/list_of_paths.txt ) | \
cat ${TMPDIR}/header.txt - > ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
${WRKDIR}/bin/rCNVmap/bin/summarize_master_burden_file.R \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2
mv ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2 \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed

#####Cut master table of bins that are significant per disease
for group in DD SCZ DD_SCZ CNCR ANY_DISEASE; do
  echo ${group}
  if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${group} ]; then
    rm -rf ${WRKDIR}/analysis/Final_Loci/significant/${group}
  fi
  mkdir ${WRKDIR}/analysis/Final_Loci/significant/${group}
  for CNV in CNV DEL DUP ANY_CNV; do
    echo ${CNV}
    for filt in all coding noncoding; do
      echo ${filt}
      list=`mktemp`
      echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
      ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
      ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      rm ${list}
      ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
      bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct \
      -i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed > \
      ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
    done
    if [ ${CNV} == "ANY_CNV" ] || [ ${group} == "ANY_DISEASE" ]; then
      filt=ANY_FILTER
      echo ${filt}
      list=`mktemp`
      echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
      ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
      ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      rm ${list}
      ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
      bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct \
      -i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed > \
      ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
    fi
  done
done
#Print table
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      echo -e "${group}\n${CNV}\n${filt}"
      #Count total number of bins tested
      zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
      #Count number of nominally significant bins from burden test
      zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | awk '{ if ($NF<=0.05) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
      #Count number of Bonferroni significant bins from burden test
      zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | awk '{ if ($NF<=0.05/8875) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
      #Count number of bins that passed permutation test
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      fgrep -v "#" | wc -l
      #Count number of merged loci that passed permutation test
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
      fgrep -v "#" | wc -l
    done | paste - - - - - - - -
  done
done
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      #Count raw overlaps with syndromic loci
      bedtools intersect -wa -u \
      -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
      -a /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed | wc -l
    done
  done
done | awk '{ print "=\""$1"/38\"" }'

#Prep for comparison of final loci within and between sets
for group in DD SCZ DD_SCZ CNCR; do
  #Make list of all loci
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
    done
  done | sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v group=${group} '{ print $0, group }' > ${TMPDIR}/${group}_master.bed
  #Make lists of coding/noncoding loci for dels + dups
  for filt in all coding noncoding; do
    for CNV in DEL DUP; do
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
    done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/${group}_${filt}.bed
  done
done
#Cat all master lists
for group in DD SCZ DD_SCZ CNCR; do
  cat ${TMPDIR}/${group}_master.bed
done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/all_master.bed
#Run comparison of final loci between del / dup within a given set
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      bedtools intersect -r -f 0.5 -c \
      -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
      -b ${TMPDIR}/${group}_${filt}.bed | awk '{ if ($NF==1) print $0 }' | wc -l
    done
  done
done
#Run comparison of loci specific to phenotype
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      awk -v group=${group} '{ if ($NF!=group) print $0 }' ${TMPDIR}/all_master.bed | \
      bedtools intersect -v -r -f 0.5 -b - \
      -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | wc -l
    done
  done
done



