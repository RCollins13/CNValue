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
mkdir ${WRKDIR}/analysis/EXON_CNV_pileups
mkdir ${WRKDIR}/analysis/EXON_CNV_burdens

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
    sed -e 's/^x/X/g' -e 's/^y/Y/g' -e 's/^MT/M/g' -e 's/^5_/5/g' -e 's/^16_/16/g' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.minSize_maxVF.maxVF.bed | \
    bedtools intersect -v -f 0.5 -a - \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
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
done | awk -v OFS="\t" '{ print $1, $2, $3, $4, $6, $8, $5, $7, $9 }'
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

#####Get pct of reference genome covered by at least one CNV
gsize=$( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | \
  awk -v OFS="\t" '{ print $1, "1", $2 }' | bedtools subtract -a - \
  -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed | \
  awk '{ sum+=($3-$2) }END{ print sum }' )
#All subjects
zcat ${WRKDIR}/data/CNV/CNV_MASTER/*CNV.GRCh37.bed.gz | fgrep -v "#" | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
-b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
  '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
#Controls only
zcat ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.CNV.GRCh37.bed.gz | fgrep -v "#" | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
-b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
  '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
#Cases only
zcat ${WRKDIR}/data/CNV/CNV_MASTER/DD.CNV.GRCh37.bed.gz \
${WRKDIR}/data/CNV/CNV_MASTER/SCZ.CNV.GRCh37.bed.gz \
${WRKDIR}/data/CNV/CNV_MASTER/CNCR.CNV.GRCh37.bed.gz | fgrep -v "#" | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
-b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
  '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'

#####Run TBRden pileup for all tissue types and all phenotypes (coding & noncoding CNVs)
#Create master mask of N-masked regions and 1Mb flanking telomeres/centromeres
grep -e 'centromere\|telomere' /data/talkowski/rlc47/src/GRCh37_heterochromatin.bed | \
awk -v OFS="\t" '{ print $1, $2-1000000, $3+1000000 }' | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' | \
cat - /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
<( grep -e 'X\|Y\|M' ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed 
#Run TBRden pileups
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed \
    -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed  \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed \
    -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed  \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_coding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 -d 1000000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.coding.bed \
    -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed  \
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

#####Run 10k CNV shift permutation tests for all comparisons
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
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 1000000 -N 100 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_all \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
      #Coding CNVs only
      bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_coding \
      -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
      #Noncoding CNVs only
      bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_noncoding \
      -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
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
      write.table(data.frame(\"perms_less_sig\"=10000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/10001),\
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
gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
${WRKDIR}/bin/rCNVmap/bin/summarize_master_burden_file.R \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2
mv ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2 \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed

#####Cut master table of bins that are significant per disease 
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${group} ]; then
    rm -rf ${WRKDIR}/analysis/Final_Loci/significant/${group}
  fi
  mkdir ${WRKDIR}/analysis/Final_Loci/significant/${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      echo ${filt}
      list=`mktemp`
      echo -e "${group}.${CNV}.${filt}.obs_p" > ${list}
      ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed \
      ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R -t 0.000001865532 \
      -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed \
      ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
      ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
      ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      rm ${list}
      #Make list of loci (±300kb merge distance)
      ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
      bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct -d 300000 \
      -i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed > \
      ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
    done
  done
done
#Print table
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    for filt in ANY_FILTER all coding noncoding; do
      if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz ]; then
        echo -e "${group}\n${CNV}\n${filt}"
        #Count total number of bins tested
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/DD_vs_CTRL/DD_vs_CTRL_CNV_all.TBRden_results.bed.gz | \
        sed '1d' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
        #Count number of nominally significant bins from burden test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of Bonferroni significant bins from burden test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of bins that passed permutation test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of merged loci that passed permutation test (as compared against master master master merged list of loci)
        bedtools intersect -u \
        -a <( zcat ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | cut -f1-3 ) \
        -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
      fi
    done | paste - - - - - - - -
  done
done
#Count overlap with syndromic CNV intervals
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    for filt in ANY_FILTER all coding noncoding; do
      #Count raw overlaps with syndromic loci
      bedtools intersect -u \
      -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
      bedtools intersect -wa -u -b - \
      -a /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed | wc -l
    done
  done
done | awk '{ print "=\""$1"/38\"" }'
#Print overlaps with syndromic CNV intervals for non-coding CNVs
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo -e "\n${group}\n"
  for CNV in ANY_CNV CNV DEL DUP; do
      bedtools intersect -u \
      -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      bedtools intersect -wa -u -b - \
      -a /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed
  done | sort -Vk1,1 -k2,2n -k3,3n | uniq
done

#####Print top-5 most significant noncoding loci per comparison
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo -e "\n${group}\n"
  while read chr start end; do
    col=$( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.noncoding.perm_signif_loci.bed.gz | \
    head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${group}.ANY_CNV.noncoding.obs_p" ) | cut -f1 )
    bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
    -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
    sort -nk${col},${col} | head -n1 | cut -f1-4,${col}
  done < <( bedtools intersect -u -wa \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
    cut -f1-3 | fgrep -v "#" ) | sort -nrk5,5
done

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

#####Get list of most significant windows per locus per comparison (drawn vs master loci)
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      #Parallelize
      bsub -q short -sla miket_sc -u nobody -J ${group}.${CNV}.${filt}.findPeakWindows \
      "${WRKDIR}/bin/rCNVmap/bin/find_peak_window.sh ${group} ${CNV} ${filt}"
      # while read chr start end; do
      #   col=$( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
      #   head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${group}.${CNV}.${filt}.obs_p" ) | cut -f1 )
      #   bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
      #   -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      #   sort -nk${col},${col} | head -n1 | cut -f1-4
      # done < <( bedtools intersect -u -wa \
      #   -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      #   -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      #   cut -f1-3 | fgrep -v "#" ) > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
    done
  done
done
#Copy to directory for plotting
if [ -e ${WRKDIR}/data/plot_data/signif_peak_windows ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_peak_windows
fi
mkdir ${WRKDIR}/data/plot_data/signif_peak_windows
for group in DD SCZ DD_SCZ CNCR ANY_DISEASE; do
  for CNV in CNV DEL DUP ANY_CNV; do
    for filt in all coding noncoding; do
      cp ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed.gz \
      ${WRKDIR}/data/plot_data/signif_peak_windows/
    done
  done
done

#####Get list of most representative windows per locus per comparison (drawn vs disease-specific loci)
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      #Parallelize
      bsub -q short -sla miket_sc -u nobody -J ${group}.${CNV}.${filt}.findRepWindows \
      "${WRKDIR}/bin/rCNVmap/bin/find_representative_windows.sh ${group} ${CNV} ${filt}"
      # while read chr start end; do
      #   size=$((${end}-${start}))
      #   mod=$( expr ${size} % 100000 )
      #   case ${mod} in
      #     0)
      #       paste <( seq ${start} 100000 $((${end}-100000)) ) \
      #       <( seq $((${start}+100000)) 100000 ${end} ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     25000)
      #       paste <( seq ${start} 100000 $((${end}-125000)) ) \
      #       <( seq $((${start}+100000)) 100000 $((${end}-25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     50000)
      #       paste <( seq $((${start}-25000)) 100000 $((${end}-75000)) ) \
      #       <( seq $((${start}+75000)) 100000 $((${end}+25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     75000)
      #       paste <( seq ${start} 100000 $((${end}-75000)) ) \
      #       <( seq $((${start}+100000)) 100000 $((${end}+25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #   esac
      # done < <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
      #   fgrep -v "#" | cut -f1-3 ) | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
    done
  done
done

#####Get counts of significant loci per comparison for Figure 1d barplot
for group in DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for filt in ANY_FILTER coding noncoding; do
    #Any CNV
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | cut -f1-3 )| wc -l
    #DEL-specific
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
    bedtools intersect -v -a - \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
    #DUP-specific
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
    bedtools intersect -v -a - \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
  done | paste - - -
done | paste - - - - > \
${WRKDIR}/data/plot_data/signif_loci_breakdown.counts.txt

#####Get overlap for 3-way venn between CNCR, DD, SCZ
for filt in ANY_FILTER coding noncoding; do
  echo -e "\n${filt}"
  for group in DD SCZ CNCR; do
    echo "${group}"
    notA=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${group} | head -n1 )
    notB=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${group} | tail -n1 )
    #Total
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Unique to group
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with A but not B
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with B but not A
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with A and B
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
  done | paste - - - - - -
done

#####Data for contrasts between deletion and duplication hotspots (Fig1e)
mkdir ${WRKDIR}/data/plot_data/del_vs_dup
#Get list of loci that are DEL only
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -v -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DEL_only.loci.bed
bedtools intersect -u -b ${TMPDIR}/DEL_only.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DEL_only.peak_windows.bed
#Get list of loci that are DUP only
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -v -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DUP_only.loci.bed
bedtools intersect -u -b ${TMPDIR}/DUP_only.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DUP_only.peak_windows.bed
#Get list of loci that are DEL and DUP
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -u -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DEL_and_DUP.loci.bed
bedtools intersect -u -b ${TMPDIR}/DEL_and_DUP.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DEL_and_DUP.peak_windows.bed
#Calculate average sizes
for CNV in DEL_only DUP_only DEL_and_DUP; do
  awk '{ print $3-$2 }' ${TMPDIR}/${CNV}.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.size.txt
done
#Fraction of DEL-only and DUP-only hotspots per disease
for group in CNCR DD_SCZ; do
  for dummy in 1; do
    echo ${group}
    #ALL
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
    #DEL only
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
    #DUP only
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
  done | paste - - - -
done > ${WRKDIR}/data/plot_data/del_vs_dup/del_vs_dup.count_by_class.txt
#Genes per 100kb per hotspot
for CNV in DEL_only DUP_only DEL_and_DUP; do
  while read chr start end skip; do
    size=$((${end}-${start}))
    fgrep -wf ${SFARI_ANNO}/genelists/Gencode_proteinCoding.genes.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
    bedtools intersect -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | \
    awk '$4 !~ /\-AS1/ { print $4 }' | sort | uniq | wc -l | \
    awk -v size=${size} '{ print (100000*$1)/size }'
  done < ${TMPDIR}/${CNV}.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.pc_genes_per_100kb.txt
done
#SD coverage per 100kb per hotspot
for CNV in DEL_only DUP_only DEL_and_DUP; do
  bedtools coverage -b <( cut -f1-3 ${TMPDIR}/${CNV}.loci.bed ) \
  -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed ) | \
  cut -f7 > ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.SD_coverage.txt
done
#Get max OR per hotspot per CNV class
for CNV in DEL DUP; do
  while read chr start end skip; do
    for group in CNCR DD SCZ; do
      for filt in all coding noncoding; do
        bedtools intersect -wa \
        -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | sed '1d' ) \
        -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
      done
    done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
  done < ${TMPDIR}/${CNV}_only.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.peakORs.txt
done
while read chr start end skip; do
  for group in CNCR DD SCZ; do
    for filt in all coding noncoding; do
      bedtools intersect -wa \
      -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_CNV_${filt}.TBRden_results.bed.gz | sed '1d' ) \
      -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
    done
  done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
done < ${TMPDIR}/DEL_and_DUP.loci.bed > \
${WRKDIR}/data/plot_data/del_vs_dup/DEL_and_DUP.peakORs.txt

#####Gather data for positive control functional enrichments (Fig1f)
#Get ORs
for dummy in 1; do
  echo "group"
  for anno in ClinGen_pathoCNVs ClinGen_HI_genes_count Petrovski2013_RVIS_1_genes_count \
    essential_genes_count DDD_2016_genes_count Sanders2015_TADAq_0.3_genes_count COSMIC_all_genes_count \
    AR_genes_count genes_protein_coding genes_all; do
      echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in ClinGen_pathoCNVs ClinGen_HI_genes_count Petrovski2013_RVIS_1_genes_count \
      essential_genes_count DDD_2016_genes_count Sanders2015_TADAq_0.3_genes_count COSMIC_all_genes_count \
      AR_genes_count genes_protein_coding genes_all; do
      tail -n1 ${WRKDIR}/analysis/Functional_Enrichments/${group}_ANY_CNV_coding/${anno}/${group}_ANY_CNV_coding_${anno}.annoBurden_results.txt | awk '{ print $1/$2 }'
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt

#####Prepare annotation files for functional enrichment analyses
#Make master bin file
zcat ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz | \
cut -f1-4 | sed '1d' > ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
#Subset master bin file to autosomes only
zcat ${WRKDIR}/data/annotations/GRCh37.master_bins.bed.gz | fgrep -v "X" | fgrep -v "Y" | sed '1d' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
#Genes (all) count
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed )  | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.genes_all.bed
#Genes (all) nearest
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL\|^M\|^X\|^Y" ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed ) | cut -f1-4,9 > \
${WRKDIR}/data/annotations/GRCh37.autosomes.genes_all.nearest.bed
#Genes (protein-coding) count
fgrep "protein_coding" ${SFARI_ANNO}/gencode/gencode.v25lift37.annotation.GRCh37.gtf | \
awk '{ if ($3=="gene") print $0 }' | sed 's/\;/\n/g' | fgrep "gene_name" | awk '{ print $NF }' | \
tr -d "\"" | sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
grep -v -e "^GL" | bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.genes_protein_coding.bed
#Genes (protein-coding) nearest
fgrep "protein_coding" ${SFARI_ANNO}/gencode/gencode.v25lift37.annotation.GRCh37.gtf | \
awk '{ if ($3=="gene") print $0 }' | sed 's/\;/\n/g' | fgrep "gene_name" | awk '{ print $NF }' | \
tr -d "\"" | sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
grep -v -e "^GL\|^M\|X\|^Y" | bedtools closest -d -t first -b - \
-a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.genes_protein_coding.nearest.bed
#Exons (all) count
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_count.bed
#Exons (all) nonredundant bases
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.merged.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exonic_bases.bed
#Exons (all) nearest
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL\|^X\|^Y\|^M" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.bed | \
  sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exons_nearest.bed
#Exons (protein-coding) count
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_proteinCoding_count.bed
#Exons (protein-coding) nonredundant bases
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exonic_proteinCoding_bases.bed
#Exons (protein-coding) nearest
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL\|^X\|^Y\|^M" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed | \
  sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_proteinCoding_nearest.bed
#GC content
bedtools nuc -fi ${h37} -bed ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
cut -f1-5 | sed '1d' > ${WRKDIR}/data/annotations/GRCh37.autosomes.GC_pct.bed
#Essential genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Wang2015_essential.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.essential_genes_distance.bed
#Essential genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Wang2015_essential.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.essential_genes_count.bed
#Constrained genes (RVIS + pLI union) count, tier 1 (10%)
cat ${SFARI_ANNO}/genelists/Petrovski2013_RVIS_10.genes.list ${SFARI_ANNO}/genelists/Lek2015_pLI_0.9.genes.list | \
sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.intolerant_tier1_genes_count.bed
#Constrained genes (RVIS + pLI union) count, tier 2 (5% RVIS, 1% pLI)
cat ${SFARI_ANNO}/genelists/Petrovski2013_RVIS_5.genes.list ${SFARI_ANNO}/genelists/Lek2015_pLI_0.99.genes.list | \
sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.intolerant_tier2_genes_count.bed
#Constrained genes (RVIS + pLI union) distance, tier 1 (10%) distance
cat ${SFARI_ANNO}/genelists/Petrovski2013_RVIS_10.genes.list ${SFARI_ANNO}/genelists/Lek2015_pLI_0.9.genes.list | \
sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.intolerant_tier1_genes_distance.bed
#Constrained genes (RVIS + pLI union) distance, tier 2 (5% RVIS, 1% pLI) distance
cat ${SFARI_ANNO}/genelists/Petrovski2013_RVIS_5.genes.list ${SFARI_ANNO}/genelists/Lek2015_pLI_0.99.genes.list | \
sort | uniq | fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.intolerant_tier2_genes_distance.bed
#Autosomal dominant genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Berg2013Bekhman2008_autosomalDominant.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.AD_genes_count.bed
#Autosomal dominant genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Berg2013Bekhman2008_autosomalDominant.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.AD_genes_distance.bed
#Autosomal recessive genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Berg2013Bekhman2008_autosomalRecessive.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.AR_genes_count.bed
#Autosomal recessive genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Berg2013Bekhman2008_autosomalRecessive.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.AR_genes_distance.bed
#ClinGen HI genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Clingen2015_Haploinsufficient.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ClinGen_HI_genes_count.bed
#ClinGen HI genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Clingen2015_Haploinsufficient.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.ClinGen_HI_genes_distance.bed
#ClinVar Disease-Associated (count)
fgrep -wf ${SFARI_ANNO}/genelists/Landrum2014_clinvar_diseaseAssociated.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ClinVar_DiseaseAssociated_genes_count.bed
#ClinVar Disease-Associated (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Landrum2014_clinvar_diseaseAssociated.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.ClinVar_DiseaseAssociated_genes_distance.bed
#GWAS Catalogue Genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Welter2014_GWAScatalog.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.GWAS_genes_count.bed
#GWAS Catalogue Genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Welter2014_GWAScatalog.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.GWAS_genes_distance.bed
#Tumor Suppressor Genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_TSC.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_TumorSuppressors_genes_count.bed
#Tumor Suppressor Genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_TSC.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_TumorSuppressors_genes_distance.bed
#Oncogenes (count)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_oncogene.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_Oncogenes_genes_count.bed
#Oncogenes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_oncogene.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_Oncogenes_genes_distance.bed
#COSMIC cancer genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_all.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_all_genes_count.bed
#COSMIC cancer genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/COSMIC_census_all.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_all_genes_distance.bed
#HPA Druggable (count)
fgrep -wf ${SFARI_ANNO}/genelists/Uhlen2015_HumanProteinAtlas_FDADruggableGenes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.Uhlen2015_HumanProteinAtlas_FDADruggableGenes_genes_count.bed
#HPA Druggable (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Uhlen2015_HumanProteinAtlas_FDADruggableGenes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.Uhlen2015_HumanProteinAtlas_FDADruggableGenes_genes_distance.bed
#T2D Genes (count)
fgrep -wf ${SFARI_ANNO}/genelists/Kim2010_T2D.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.Kim2010_T2D_genes_count.bed
#T2D Genes (distance)
fgrep -wf ${SFARI_ANNO}/genelists/Kim2010_T2D.genes.list \
${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Kim2010_T2D_genes_distance.bed
#Bulk SFARI gene annotations
while read list; do
  echo ${list}
  #Count
  fgrep -wf ${SFARI_ANNO}/genelists/${list}.genes.list \
  ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.${list}_genes_count.bed
  #Distance
  fgrep -wf ${SFARI_ANNO}/genelists/${list}.genes.list \
  ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
  bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
  sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${list}_genes_distance.bed
done < <( echo -e "Sugathan2014_CHD8\nPetrovski2013_RVIS_10\nPetrovski2013_RVIS_5\nPetrovski2013_RVIS_1\n\
Lek2015_pLI_0.9\nLek2015_pLI_0.99\nLek2015_pLI_0.99999\nLee2016_RBFOX1\nLiu2014_DAWN\nSanders2015_TADAq_0.1\n\
Sanders2015_TADAq_0.3\nCotney2015_CHD8_hNSC\nCotney2015_CHD8_Brain\nCHD8_union\nCHD8_intersection\n\
Gencode_FMRP\nPawson2014_GPCR\nWelter2014_GWAScatalog\nRVIS1_pLI0.99999.Constrained_HC_union\n\
RVIS10_pLI0.9.Constrained_union\nDDD_2016\nDAWN_TADAq0.3_ASD_ID_EPI_MC_DDD2016.NDD_union" )
#Old SFARI ASD, EPI & ID gene lists
for group in ASD EPI ID; do
  for level in HC MC LC; do
    #Count
    fgrep -wf /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/genes/${group}_${level}.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL" | \
    bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
    sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.${group}_${level}_genes_count.bed
    #Distance
    fgrep -wf /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/genes/${group}_${level}.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | grep -v -e "^GL\|^X\|^Y\|^M" | \
    bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -b - | \
    sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
    ${WRKDIR}/data/annotations/GRCh37.autosomes.${group}_${level}_genes_distance.bed
  done
done
#ENCODE DNAseH1 - Any
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_raw.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_raw_coverage.bed
#ENCODE DNAseH1 - 50% Cell Types
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_50pct_celltypes.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_50pct_coverage.bed
#ENCODE DNAseH1 - 90% Cell Types
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_90pct_celltypes.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_90pct_coverage.bed
#UCNEs (distance)
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -ve "^X\|^Y\|^M\|^GL"  ${SFARI_ANNO}/noncoding/Dimitrieva2013_UCNE.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > \
${WRKDIR}/data/annotations/GRCh37.autosomes.UCNE_distance.bed
#UCNEs (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/noncoding/Dimitrieva2013_UCNE.bed | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.UCNE_count.bed
#RepeatMasker by class
for class in LINE SINE LTR Satellite Simple_repeat Low_complexity; do
  bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -a <( grep -v -e "^GL" /data/talkowski/rlc47/src/GRCh37.${class}.RMSK.bed ) | \
  awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.RM_${class}_coverage.bed
done
#HARs (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/noncoding/Doan2016_HAR.bed | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.HAR_count.bed
#HARs (distance)
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -ve "^X\|^Y\|^M\|^GL" ${SFARI_ANNO}/noncoding/Doan2016_HAR.bed | sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > \
${WRKDIR}/data/annotations/GRCh37.autosomes.HAR_distance.bed
#SegDups
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.SegDup_coverage.bed
#GM12878 histone marks
for mark in H3k27ac H3k4me1 H3k4me2 H3k4me3 H2az H3k36me3 H3k9ac H3k9me3; do
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeBroadHistoneGm12878${mark}StdSig.bigWig \
  <( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab
  paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab ) > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.bed
  rm ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab
done
#ENCODE 100mer alignability & 35bp uniqueness
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeCrgMapabilityAlign100mer.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeDukeMapabilityUniqueness35bp.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab
#GM12878 nucleosome map
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeSydhNsomeGm12878Sig.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab
#ENCODE TF ChIP peaks (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | sed 's/^chr//g' | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.TF_ChIP_peaks.bed
#ENCODE TF ChIP peaks (count by TF)
mkdir ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF
zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | sed 's/^chr//g' > \
${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.all.bed
zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
awk '{ if ($5<250) print $0 }' | sed 's/^chr//g' > \
${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.0_250.bed
zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
awk '{ if ($5>=250 && $5<500) print $0 }' | sed 's/^chr//g' > \
${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.250_500.bed
zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
awk '{ if ($5>=500 && $5<750) print $0 }' | sed 's/^chr//g' > \
${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.500_750.bed
zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
awk '{ if ($5>=750) print $0 }' | sed 's/^chr//g' > \
${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.750_1000.bed
gzip ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/*
while read TF; do
  echo ${TF}
  mkdir ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}
  zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v TF=${TF} '{ if ($4==TF) print $0 }' | sed 's/^chr//g' > \
  ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.all.bed
  zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v TF=${TF} '{ if ($4==TF && $5<250) print $0 }' | sed 's/^chr//g' > \
  ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.0_250.bed
  zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v TF=${TF} '{ if ($4==TF && $5>=250 && $5<500) print $0 }' | sed 's/^chr//g' > \
  ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.250_500.bed
  zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v TF=${TF} '{ if ($4==TF && $5>=500 && $5<750) print $0 }' | sed 's/^chr//g' > \
  ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.500_750.bed
  zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | \
  awk -v TF=${TF} '{ if ($4==TF && $5>=750) print $0 }' | sed 's/^chr//g' > \
  ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.750_1000.bed
  gzip ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/*
done < <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | cut -f4 | sort | uniq )
while read TF; do
  echo ${TF}
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.all.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ENCODE_${TF}_peaks_all.bed
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.0_250.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ENCODE_${TF}_peaks_0_250.bed
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.250_500.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ENCODE_${TF}_peaks_250_500.bed
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.500_750.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ENCODE_${TF}_peaks_500_750.bed
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.750_1000.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ENCODE_${TF}_peaks_750_1000.bed
done < <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | cut -f4 | sort | uniq | \
  cat <( echo -e "all_TF" ) - )
#ClinGen Pathogenic CNV Regions (distance)
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( zcat ${SFARI_ANNO}/misc/ClinGen_PathogenicCNVs.bed.gz | grep -ve "^X\|^Y\|^M\|^GL" | sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > \
${WRKDIR}/data/annotations/GRCh37.autosomes.ClinGen_pathoCNVs_distance.bed
#ClinGen Pathogenic CNV Regions (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/misc/ClinGen_PathogenicCNVs.bed.gz | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.ClinGen_pathoCNVs.bed
#GWAS catalog variants (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/misc/GWAS_Catalog_variants.merged_simple.bed.gz | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.GWAS_loci.bed
#GTEx eQTLs (count & distance)
while read tissue; do
  echo ${tissue}
  mkdir ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}
  sed '1d' ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}_Analysis.snpgenes | \
  awk -v OFS="\t" '{ print $14, $15, $15+1 }' > \
  ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' )
while read tissue; do
  echo ${tissue}
  #Count
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed | \
  sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_eQTLs_${tissue}_count.bed
  #Distance
  bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b <( grep -ve "^X\|^Y\|^M\|^GL" ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed | \
    sort -Vk1,1 -k2,2n -k3,3n ) | sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,8 > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_eQTLs_${tissue}_distance.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' )
#COSMIC variants (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/misc/COSMIC_variants.bed.gz | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.COSMIC_variants.bed
#Affy6 Probes (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/misc/Affy6_Probes.bed.gz | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Affy6_Probes.bed
#Illumina 1MDuo Probes (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/misc/Illumina_1MDuo_Probes.bed.gz | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Illumina_1MDuo_Probes.bed
#PhyloP 100-way conservation
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
/scratch/miket/rlc47temp/DOSAGE/hg19.100way.phyloP100way.bw \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
#dbSNP 142 variants (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b /scratch/miket/rlc47temp/DOSAGE/snp142_SNVs_MNPs.bed | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.dbSNP142_vars.bed
#Replication timing
bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o mean -c 4 \
-b <( fgrep -v "X" /scratch/miket/rlc47temp/DOSAGE/hg19_repTiming.bed | fgrep -v "Y" ) | \
awk -v OFS="\t" '{ if ($4==".") $4="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.RepTiming.bed
#Recombination frequency (deCODE)
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/SexAveraged.bw \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab
#GTEx Expression
#Curate
zcat ${SFARI_ANNO}/misc/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | \
sed -n '3p' | sed 's/\t/\n/g' | sed -e 's/\-//g' -e 's/(//g' -e 's/)//g' | sed 's/ \+/_/g' | \
sed '1,2d' | paste -s | paste <( echo -e "chr\tstart\tend\tgene" ) - > \
${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed
zcat ${SFARI_ANNO}/misc/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | \
sed '1,3d' | cut -f2- | sort -k1,1 | join -t $'\t' -1 1 -2 1 \
<( awk -v OFS="\t" '{ print $4, $1, $2, $3 }' ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | sort -k1,1 ) - | \
awk -v OFS="\t" '{ if ($1!="7SK" && $1!="5S_rRNA" && $1!="Y_RNA") print $2, $3, $4, $1, $0 }' | \
cut -f 1-4,9- | sort -Vk1,1 -k2,2n -k3,3n >> \
${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed
while read tissue; do #master matrix of all genes
  echo ${tissue}
  col=$( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | fgrep -w ${tissue} | cut -f2 )
  cut -f1-4,${col} ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  sed '1d' > ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed
  gzip -f ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
while read tissue; do #Genes with FPKM > 50
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | \
  awk '{ if ($5>=50) print $0 }' > \
  ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed
  gzip -f ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
while read tissue; do #curate binwise expression levels
  echo ${tissue}
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o mean -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_avgExpression_${tissue}.bed
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o max -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_maxExpression_${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
while read tissue; do #curate counts of highly expressed genes
  echo ${tissue}
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed.gz | \
  sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_highExpressors_${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
#Schmitt compartment - Fibroblast & LCL
for cell in BL AO AD tro SX SB RV PO PA OV npc msc mes LV LI LG imr90 HC h1 GM12878 CO; do
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${cell}.pc.bw \
  <( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab
  paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab ) > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.bed
  rm ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab
done
#Schmitt TBRs
while read tissue; do
  echo -e "${tissue}"
  #Count
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${tissue}_TBRs.bed
  #Distance
  bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b <( grep -ve "^X\|^Y\|^M\|^GL" /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed | \
    sort -Vk1,1 -k2,2n -k3,3n ) | cut -f1-4,9 > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${tissue}_TBRs_distance.bed
done < <( cat <( echo "MERGED" ) /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list )
while read tissue; do
  cat /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools intersect -c -b - \
-a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.ALL_TBRs.bed
while read tissue; do
  cat /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | grep -ve "^X\|^Y\|^M\|^GL" | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools closest -d -t first -b - \
-a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | sort -Vk1,1 -k2,2n -k3,3n | \
cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.ALL_TBRs_distance.bed
#Super enhancers
while read tissue; do
  echo ${tissue}
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/superEnhancers/cleaned/${tissue}_superEnhancer.bed | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${tissue}_superEnhancer.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/superEnhancers/*bed | \
  awk '{ print $9 }' | sed 's/superEnhancers\//\t/g' | sed 's/\.bed/\t/g' | cut -f2 )
#Tissue-specific enhancers
while read tissue; do
  echo ${tissue}
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b <( cat ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed | uniq ) | \
  sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.${tissue}_enhancerAtlas.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 )
#Brain enhancer RNAs (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/noncoding/Yao2015_BrainEnhancerRNA.bed | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.BrainEnhancerRNAs.bed
#FMRP targets
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/noncoding/Ascano2012_FMRP.bed | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.FMRP_targets.bed
#FANTOM5 NPC, fetal brain, and adult brain enhancers
for tissue in NPC AdultBrain FetalBrain; do
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/FANTOM52015_${tissue}_TPM5.bed | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.FANTOM5_${tissue}_enhancer.bed
done
#Cotney CHD8 peaks
for tissue in Brain hNSC; do
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/noncoding/Cotney2015_CHD8_${tissue}_peaks.bed | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.CHD8_${tissue}_targets.bed
done
#Talkowski 36 rCNV loci (count)
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.HCrCNVs.bed
#Talkowski 36 rCNV loci (distance)
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -ve "^X\|^Y" /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.HCrCNVs_distance.bed

#####Check number of columns for each annotation file
while read anno skip; do
  zcat ${WRKDIR}/data/annotations/GRCh37.autosomes.${anno}.bed.gz | \
  head -n1 | awk '{ print NF }'
  echo "${anno}"
done < ${WRKDIR}/lists/annotations.all.list | paste - - | sort -nk1,1

#####Launch functional enrichment analyses of all significant windows
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    for filt in ANY_FILTER all coding noncoding; do
      if [ -e ${WRKDIR}/analysis/Functional_Enrichments/${group}_${CNV}_${filt} ]; then
        rm -r ${WRKDIR}/analysis/Functional_Enrichments/${group}_${CNV}_${filt}
      fi
      mkdir ${WRKDIR}/analysis/Functional_Enrichments/${group}_${CNV}_${filt}
      while read anno measure; do
        mkdir ${WRKDIR}/analysis/Functional_Enrichments/${group}_${CNV}_${filt}/${anno}
        bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_${filt} \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_annoBurden.R \
        ${WRKDIR}/data/annotations/GRCh37.autosomes.${anno}.bed.gz \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed.gz \
        ${measure} ${WRKDIR}/analysis/Functional_Enrichments/${group}_${CNV}_${filt}/${anno} \
        ${group}_${CNV}_${filt}_${anno}"
      done < ${WRKDIR}/lists/annotations.all.list
    done
  done
done

#####Recreate functional annotation tracks after restricting out elements near protein-coding exons
#Make merged list of protein-coding exons
fgrep -wf ${SFARI_ANNO}/genelists/Gencode_proteinCoding.genes.list \
${SFARI_ANNO}/gencode/gencode.v25lift37.exons_minus_UTRs.lex.bed | \
bedtools merge -c 4 -o distinct -i - | sort -Vk1,1 -k2,2n -k3,3n > \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed
#Make list of protein-coding exons with ±10kb flanks
awk -v OFS="\t" '{ print $1, $2-10000, $3+10000 }' \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.merged.bed | \
awk -v OFS="\t" '{ if ($2<1) $2=1; print }' | \
awk -v OFS="\t" '{ if ($3<2) $3=2; print }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed
#Make list of subtracted master bins
bedtools subtract -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
awk -v OFS="\t" '{ print "chr"$0, NR }' > ${TMPDIR}/Subtracted_bins.bed
#GC content
bedtools subtract -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools nuc -fi ${h37} -bed - | cut -f1-5 | sed '1d' > ${TMPDIR}/GC.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/GC.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.GC_pct.bed
#ENCODE DNAseH1 - Any
bedtools intersect -v -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_raw.bed ) \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a - | awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.DNAseH1_raw_coverage.bed
#ENCODE DNAseH1 - 50% Cell Types
bedtools intersect -v -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_50pct_celltypes.bed ) \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a - | awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.DNAseH1_50pct_coverage.bed
#ENCODE DNAseH1 - 90% Cell Types
bedtools intersect -v -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_90pct_celltypes.bed ) \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a - | awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.DNAseH1_90pct_coverage.bed
#UCNEs
bedtools intersect -v -a ${SFARI_ANNO}/noncoding/Dimitrieva2013_UCNE.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.UCNE_count.bed
#HARs
bedtools intersect -v -a ${SFARI_ANNO}/noncoding/Doan2016_HAR.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.HAR_count.bed
#RepeatMasker by class
for class in LINE SINE LTR Satellite Simple_repeat Low_complexity; do
  echo "${class}"
  bedtools intersect -v -a <( grep -v -e "^GL" /data/talkowski/rlc47/src/GRCh37.${class}.RMSK.bed ) \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -a - | awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.RM_${class}_coverage.bed
done
#SegDups
bedtools intersect -v -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed ) \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a - | awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.SegDup_coverage.bed
#GM12878 histone marks
for mark in H3k27ac H3k4me1 H3k4me2 H3k4me3 H2az H3k36me3 H3k9ac H3k9me3; do
  echo "${mark}"
  bedtools subtract -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  awk -v OFS="\t" '{ print "chr"$0, NR }' > ${TMPDIR}/${mark}.originalbins.bed
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeBroadHistoneGm12878${mark}StdSig.bigWig \
  <( cut -f1-3,5 ${TMPDIR}/${mark}.originalbins.bed ) \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${mark}.tab
  paste <( cut -f1-4 ${TMPDIR}/${mark}.originalbins.bed | sed 's/chr//g' ) \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${mark}.tab ) > \
  ${TMPDIR}/${mark}.subtracted.bed
  ${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
  ${TMPDIR}/${mark}.subtracted.bed \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${mark}.bed
  rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${mark}.tab
done
#ENCODE 100mer alignability & 35bp uniqueness
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeCrgMapabilityAlign100mer.bigWig \
<( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Align_100mer.tab
paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Align_100mer.tab ) > \
${TMPDIR}/Align_100mer.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/Align_100mer.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Align_100mer.bed
rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Align_100mer.tab
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeDukeMapabilityUniqueness35bp.bigWig \
<( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Unique_35bp.tab
paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Unique_35bp.tab ) > \
${TMPDIR}/Unique_35bp.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/Unique_35bp.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Unique_35bp.bed
rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Unique_35bp.tab
#GM12878 nucleosome map
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeSydhNsomeGm12878Sig.bigWig \
<( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Nucleosome_Density.tab
paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Nucleosome_Density.tab ) > \
${TMPDIR}/Nucleosome_Density.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/Nucleosome_Density.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Nucleosome_Density.bed
rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Nucleosome_Density.tab
#ENCODE TF ChIP peaks (count)
bedtools intersect -v -a <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | sed 's/^chr//g' | cut -f1-3 ) \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.TF_ChIP_peaks.bed
#ENCODE TF ChIP peaks (count by TF)
while read TF; do
  echo ${TF}
  for slice in all 0_250 250_500 500_750 750_1000; do
    bedtools intersect -v -a ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.${slice}.bed.gz \
    -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
    bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
    -b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.ENCODE_${TF}_peaks_${slice}.bed
  done
done < <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | cut -f4 | sort | uniq | \
  cat <( echo -e "all_TF" ) - )
#GWAS catalog variants (count)
bedtools intersect -v -a ${SFARI_ANNO}/misc/GWAS_Catalog_variants.merged_simple.bed.gz \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.GWAS_loci.bed
#GTEx eQTLs (count & distance)
while read tissue; do
  echo ${tissue}
  #Count
  bedtools intersect -v -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.GTEx_eQTLs_${tissue}_count.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' )
#COSMIC variants (count)
bedtools intersect -v -a ${SFARI_ANNO}/misc/COSMIC_variants.bed.gz \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.COSMIC_variants.bed
#Affy6 Probes (count)
bedtools intersect -v -a ${SFARI_ANNO}/misc/Affy6_Probes.bed.gz \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Affy6_Probes.bed
#Illumina 1MDuo Probes (count)
bedtools intersect -v -a ${SFARI_ANNO}/misc/Illumina_1MDuo_Probes.bed.gz \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Illumina_1MDuo_Probes.bed
#PhyloP 100-way conservation
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
/scratch/miket/rlc47temp/DOSAGE/hg19.100way.phyloP100way.bw \
<( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.PhyloP_Conservation.tab
paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.PhyloP_Conservation.tab ) > \
${TMPDIR}/PhyloP_Conservation.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/PhyloP_Conservation.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.PhyloP_Conservation.bed
rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.PhyloP_Conservation.tab
#dbSNP 142 variants (count)
bedtools intersect -v -a /scratch/miket/rlc47temp/DOSAGE/snp142_SNVs_MNPs.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.dbSNP142_vars.bed
#Replication timing
bedtools map -a <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' | sort -Vk1,1 -k2,2n -k3,3n ) -o mean -c 4 \
-b <( fgrep -v "X" /scratch/miket/rlc47temp/DOSAGE/hg19_repTiming.bed | fgrep -v "Y" ) | \
awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' > ${TMPDIR}/hg19_repTiming.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/hg19_repTiming.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.RepTiming.bed
#Recombination frequency (deCODE)
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/SexAveraged.bw \
<( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Recomb_Freq.tab
paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Recomb_Freq.tab ) > \
${TMPDIR}/Recomb_Freq.subtracted.bed
${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
${TMPDIR}/Recomb_Freq.subtracted.bed \
${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Recomb_Freq.bed
rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.Recomb_Freq.tab
#GTEx Expression
while read tissue; do
  echo ${tissue}
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o mean -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.GTEx_avgExpression_${tissue}.bed
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o max -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.GTEx_maxExpression_${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
#Schmitt compartment - Fibroblast & LCL
for cell in BL AO AD tro SX SB RV PO PA OV npc msc mes LV LI LG imr90 HC h1 GM12878 CO; do
  echo ${cell}
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${cell}.pc.bw \
  <( cut -f1-3,5 ${TMPDIR}/Subtracted_bins.bed ) \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${cell}_compartment.tab
  paste <( cut -f1-4 ${TMPDIR}/Subtracted_bins.bed | sed 's/chr//g' ) \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${cell}_compartment.tab ) > \
  ${TMPDIR}/${cell}_compartment.subtracted.bed
  ${WRKDIR}/bin/rCNVmap/bin/subtractedBed_weightedMeans.R \
  ${TMPDIR}/${cell}_compartment.subtracted.bed \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${cell}_compartment.bed
  rm ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${cell}_compartment.tab
done
#Schmitt TBRs
while read tissue; do
  #Count
  bedtools intersect -v -a /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${tissue}_TBRs.bed
  #Distance
  grep -ve "^X\|^Y\|^M\|^GL" /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed | \
  bedtools intersect -v -a - -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools closest -d -t first -b - -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | cut -f1-4,9 > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${tissue}_TBRs_distance.bed
done < <( cat <( echo "MERGED" ) /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list )
while read tissue; do
  cat /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools intersect -v -a - \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -b - \
-a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.ALL_TBRs.bed
while read tissue; do
  cat /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | grep -ve "^X\|^Y\|^M\|^GL" | \
bedtools intersect -v -a - \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools closest -d -t first -b - \
-a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | sort -Vk1,1 -k2,2n -k3,3n | \
cut -f1-4,9 > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.ALL_TBRs_distance.bed
#Super enhancers
while read tissue; do
  echo ${tissue}
  #Count
  bedtools intersect -v -a ${SFARI_ANNO}/noncoding/superEnhancers/cleaned/${tissue}_superEnhancer.bed \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${tissue}_superEnhancer_count.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/superEnhancers/*bed | \
  awk '{ print $9 }' | sed 's/superEnhancers\//\t/g' | sed 's/\.bed/\t/g' | cut -f2 )
#Tissue-specific enhancers
while read tissue; do
  echo ${tissue}
  cat ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed | sort -Vk1,1 -k2,2n -k3,3n | uniq | \
  bedtools intersect -v -a - -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${tissue}_enhancerAtlas.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 )
#All enhancers
cat ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/*_enhancers.bed | \
bedtools intersect -v -a - -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.ALL_TISSUES_enhancerAtlas.bed
#High-confidence enhancers (overlap with a strong TFBS)
while read tissue; do
  echo ${tissue}
  cat ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed | sort -Vk1,1 -k2,2n -k3,3n | uniq | \
  bedtools intersect -v -a - -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -u -a - -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.750_1000.bed.gz | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${tissue}_enhancerAtlas_strong.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 )
cat ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/*_enhancers.bed | \
bedtools intersect -v -a - -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -u -a - -b ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/all_TF/all_TF.ENCODE_peaks.750_1000.bed.gz | \
sort -Vk1,1 -k2,2n -k3,3n | uniq | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.ALL_TISSUES_enhancerAtlas_strong.bed
#Brain enhancer RNAs
bedtools intersect -v -a ${SFARI_ANNO}/noncoding/Yao2015_BrainEnhancerRNA.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.BrainEnhancerRNAs.bed
#FMRP targets
bedtools intersect -v -a ${SFARI_ANNO}/noncoding/Ascano2012_FMRP.bed \
-b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b - | sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.FMRP_targets.bed
#FANTOM5 NPC, fetal brain, and adult brain enhancers
for tissue in NPC AdultBrain FetalBrain; do
  bedtools intersect -v -a ${SFARI_ANNO}/noncoding/FANTOM52015_${tissue}_TPM5.bed \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.FANTOM5_${tissue}_enhancer.bed
done
#Cotney CHD8 peaks
for tissue in Brain hNSC; do
  bedtools intersect -v -a ${SFARI_ANNO}/noncoding/Cotney2015_CHD8_${tissue}_peaks.bed \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed | \
  bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
  -b - | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.CHD8_${tissue}_targets.bed
done

#####Launch exon-excluded functional enrichment analyses of noncoding significant windows
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_${CNV}_noncoding ]; then
      rm -r ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_${CNV}_noncoding
    fi
    mkdir ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_${CNV}_noncoding
    while read anno measure; do
      mkdir ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_${CNV}_noncoding/${anno}
      bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_noncoding_exonExclusion \
      "${WRKDIR}/bin/rCNVmap/bin/TBRden_annoBurden.R \
      ${WRKDIR}/data/annotations/exonExclusion/GRCh37.autosomes.noExons.${anno}.bed.gz \
      ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.noncoding.perm_signif_bins.bed.gz \
      ${measure} ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_${CNV}_noncoding/${anno} \
      ${group}_${CNV}_noncoding_${anno}"
    done < ${WRKDIR}/lists/annotations_exonExclusion.all.list
  done
done

#####Collect master list of noncoding enrichments
while read anno skip; do
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      if [ -e ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_ANY_CNV_noncoding/${anno}/${group}_ANY_CNV_noncoding_${anno}.annoBurden_results.txt ]; then
        tail -n1 ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion/${group}_ANY_CNV_noncoding/${anno}/${group}_ANY_CNV_noncoding_${anno}.annoBurden_results.txt | \
        awk -v OFS="\t" '{ print $1/$2, $4 }' 2> /dev/null
      else
        echo -e "NA\nNA"
      fi
    done | paste -s 
  done | paste -s
done < ${WRKDIR}/lists/annotations_exonExclusion.all.list > \
${TMPDIR}/noncoding_anno_results.txt

#####Run protein-coding exon pileups for noncoding rCNVs in flanks
#Prep list of 50kb flanks
distance=50000
awk -v OFS="\t" -v d=${distance} '{ print $1, $2-d, $3+d, $4 }' \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed > \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.50kb_windows.bed
#Launch tests
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_exon_flanking_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 10000 -n -z -t upper \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    -p ${group}_${CNV}_noncoding_exonFlanks \
    -o ${WRKDIR}/analysis/EXON_CNV_burdens/${group}.${CNV}.TBRden_exon_burdens.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.noncoding.bed.gz
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.50kb_windows.bed \
    ${WRKDIR}/analysis/EXON_CNV_burdens/"
  done
done

#####Run burden test on protein-coding exon pileups
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_exon_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/EXON_CNV_pileups/CTRL.${CNV}.TBRden_exon_pileup.bed.gz \
    ${WRKDIR}/analysis/EXON_CNV_pileups/${group}.${CNV}.TBRden_exon_pileup.bed.gz \
    ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_exons ${color}"
  done
done

# #####Run 100k CNV shift permutation tests for exon comparisons
# #Note: initial p-value cutoff used: 0.05/18991 = 0.000002632826
# #This corresponds to the number of non-overlapping protein-coding exons we tested
# for group in DD SCZ DD_SCZ CNCR; do
#   echo ${group}
#   if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL ]; then
#     rm -rf ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
#   fi
#   mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
#   mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split
#   for CNV in DEL DUP CNV; do
#     echo ${CNV}
#     #Coding + noncoding CNVs
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
#     #Coding CNVs only
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed
#     #Noncoding CNVs only
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
#     #Split into partitions of 1k permutations each (x100 per comparison)
#     for i in $( seq -w 001 100 ); do
#       echo ${i}
#       #Coding + noncoding CNVs
#       if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i} ]; then
#         mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}
#       fi
#       bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
#       "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 1000000 -N 1000 \
#       -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
#       -p ${group}_vs_CTRL_${CNV}_all \
#       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
#       #Coding CNVs only
#       bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
#       "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 1000 \
#       -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
#       -p ${group}_vs_CTRL_${CNV}_coding \
#       -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
#       #Noncoding CNVs only
#       bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
#       "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 1000 \
#       -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
#       -p ${group}_vs_CTRL_${CNV}_noncoding \
#       -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#       ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
#     done
#   done
# done






