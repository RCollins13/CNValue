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
mkdir ${WRKDIR}/bin/LSF
mkdir ${WRKDIR}/data
mkdir ${WRKDIR}/data/CNV
mkdir ${WRKDIR}/data/CNV/CNV_RAW
mkdir ${WRKDIR}/data/CNV/CNV_MASTER
mkdir ${WRKDIR}/data/CNV/CNV_DENSITY
mkdir ${WRKDIR}/data/plot_data
mkdir ${WRKDIR}/data/unfiltered_annotations
mkdir ${WRKDIR}/data/filtered_annotations
mkdir ${WRKDIR}/data/misc
mkdir ${WRKDIR}/analysis
mkdir ${WRKDIR}/analysis/CNV_set_enrichments
mkdir ${WRKDIR}/analysis/BIN_CNV_pileups
mkdir ${WRKDIR}/analysis/BIN_CNV_burdens
mkdir ${WRKDIR}/analysis/BIN_CNV_permutation
mkdir ${WRKDIR}/analysis/Final_Loci
mkdir ${WRKDIR}/analysis/Final_Loci/significant
mkdir ${WRKDIR}/analysis/hotspot_enrichments
mkdir ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations
mkdir ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations
mkdir ${WRKDIR}/analysis/EXON_CNV_burdens
mkdir ${WRKDIR}/analysis/TBR_CNV_pileups
mkdir ${WRKDIR}/analysis/TBR_CNV_burdens

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
while read group tier description HPO_in HPO_out; do
  echo -e "${group}"
  echo -e "${HPO_in}" | sed 's/\;/\n/g' | fgrep -wf - \
  <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
  fgrep -wvf <( echo -e "${HPO_out}" | sed 's/\;/\n/g' ) | wc -l
  echo -e "${description}"
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | paste - - -

#####Assign HPO terms to all CNVs
#Cut HPO mappings 
cut -f1,3 ${WRKDIR}/data/HPO_map/HPO_map.modified.tsv > \
${TMPDIR}/HPO_mappings_for_application.txt
#Create list of HPO mappings for all germline CNVs
for cohort in SSC Coe Talkowski PGC; do
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
for cohort in SSC Coe Talkowski PGC; do
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
for cohort in SSC Coe Talkowski PGC; do
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
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV
for CNV in DEL DUP; do
  #Create master list of all CNVs
  for cohort in Shaikh Suktitipat Uddin Vogler Cooper SSC Coe Talkowski PGC TCGA; do
    if [ -e ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz ]; then
      zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_GERM.${CNV}.raw.bed.gz | fgrep -v "#" | \
      awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
    fi
    if [ -e ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CTRL.${CNV}.raw.bed.gz ]; then
      zcat ${WRKDIR}/data/CNV/CNV_RAW/${cohort}_CNVs/${cohort}_CTRL.${CNV}.raw.bed.gz | fgrep -v "#" | \
      awk -v OFS="\t" '{ print $1, $2, $3, $4, $4, $5 }'
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

#####Filter merged CNVs across cohort (max VF)
if [ -e ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV ]; then
  rm -rf ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/
fi
mkdir ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/
max_VF=0.01
for CNV in DEL DUP; do
  #Germline
  awk -v max_VF=${max_VF} '{ if (($9/122312)<max_VF) print $0 }' \
  ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/germline_${CNV}.merged.maxVF.bed
  #Cancer
  cat ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/cancer_${CNV}.merged.bed > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/cancer_${CNV}.merged.maxVF.bed
done

#####Filter merged CNVs per study (max VF per study) & reassign original coordinates
max_VF=0.01
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
    study_base=$( echo "${study}" | cut -f1 -d_ )
    pheno=$( echo ${study} | sed 's/_/\t/g' | awk '{ print $NF }' )
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*.maxVF.bed | fgrep "${study}_${CNV}" | \
    cut -f7 | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
    awk -v max_VF=${max_VF} -v n=${n} '{ if (($2/n)<max_VF) print $1 }' | \
    fgrep -wf - ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*.maxVF.bed | \
    cut -f2- -d\: | cut -f4 | fgrep "${study}_${CNV}" | sort | uniq | \
    fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/${study_base}_CNVs/${study}.${CNV}.raw.bed.gz ) | \
    awk -v OFS="\t" -v PMID=${PMID} '{ print $1, $2, $3, $4, $5, $6, PMID }' | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  done < <( fgrep -v "TCGA_CNCR" ${WRKDIR}/lists/Studies_SampleSizes.list )
  echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
  fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/cancer_${CNV}.merged.maxVF.bed | \
  cut -f4 | fgrep -wf - <( zcat ${WRKDIR}/data/CNV/CNV_RAW/TCGA_CNVs/TCGA_CNCR.${CNV}.raw.bed.gz ) >> \
  ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/TCGA_CNCR_${CNV}.merged.maxVF.maxVF.originalCoords.bed
done

#####Filter merged CNVs on minimum size
min_size=20000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
    awk -v min_size=${min_size} '{ if ($3-$2>min_size) print $0 }' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.bed | \
    bedtools intersect -v -f 0.5 -a - -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Filter merged CNVs on exclusion loci
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
    sed -e 's/^x/X/g' -e 's/^y/Y/g' -e 's/^MT/M/g' -e 's/^5_/5/g' -e 's/^16_/16/g' \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.bed | \
    bedtools intersect -v -f 0.5 -a - \
    -b ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Filter merged CNVs on maximum size
max_size=5000000
for CNV in DEL DUP; do
  while read study n PMID; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
    fgrep -v "#" ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
    awk -v max_size=${max_size} '{ if ($3-$2<max_size) print $0 }' >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/${study}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed
  done < ${WRKDIR}/lists/Studies_SampleSizes.list
done

#####Make master germline, control, and cancer files (post-filter)
for pheno in GERM CTRL CNCR; do
  for CNV in DEL DUP; do
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${pheno}_${CNV}.final_filters.bed
    echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.maxSize.bed | \
    fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${pheno}_${CNV}.final_filters.bed
    cat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*${pheno}_${CNV}.merged.maxVF.maxVF.originalCoords.minSize.blacklisted.bed | \
    fgrep -v "#" | sort -Vk1,1 -k2,2n -k3,3n >> \
    ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${pheno}_${CNV}.final_filters_noMaxSize.bed
  done
done
#Gzip all filtered CNV files
gzip ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/*

#####Final merger per analysis group, filtered on max size & unfiltered on max size
if [ -e ${WRKDIR}/data/CNV/CNV_MASTER/ ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_MASTER/
fi
mkdir ${WRKDIR}/data/CNV/CNV_MASTER/
while read group eti tier descrip include exclude color; do
  if [ ${exclude} != "NA" ]; then
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -wvf <( echo ${exclude} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
    done
  else
    for CNV in DEL DUP; do
      #Filtered on max size (main CNV set)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${eti}_${CNV}.final_filters.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed
      #Unfiltered on max size (used for size distribs)
      echo -e "#chr\tstart\tend\tVID\tCNV\tPheno\tSource_PMID" > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
      zcat ${WRKDIR}/data/CNV/CNV_RAW/filtered_CNV/ALL_${eti}_${CNV}.final_filters_noMaxSize.bed.gz | \
      fgrep -wf <( echo ${include} | sed 's/\;/\n/g' ) | \
      fgrep -v "#" | sed -e 's/^5_/5/g' -e 's/^y/Y/g' -e 's/^16_/16/g' -e 's/^x_/X/g' -e 's/^23/X/g' -e 's/^24/Y/g' | \
      sort -Vk1,1 -k2,2n -k3,3n >> ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed
    done
  fi
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Make merged CNV (DEL+DUP) set
while read group eti tier descrip include exclude color; do
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
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Generate coding and noncoding CNV callsets
while read group eti tier descrip include exclude color; do
  echo ${group}
  for class in CNV DEL DUP; do
    echo ${class}
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
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )
gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/*bed

#####Generate CNV density per callset
if [ -e ${WRKDIR}/data/CNV/CNV_DENSITY ]; then
  rm -r ${WRKDIR}/data/CNV/CNV_DENSITY
fi
mkdir ${WRKDIR}/data/CNV/CNV_DENSITY
while read group eti tier descrip include exclude color; do
  echo ${group}
  bsub -q normal -sla miket_sc -J ${group}_CNV_densities -u nobody \
  "${WRKDIR}/bin/rCNVmap/bin/generate_CNV_densities.sh ${group}"
  # for CNV in DEL DUP CNV; do
  #   echo ${CNV}
  #   bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  #   -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz > \
  #   ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
  #   gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.all.CNV_density.bg
  #   bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  #   -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz > \
  #   ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
  #   gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.noncoding.CNV_density.bg
  #   bedtools genomecov -bga -g /data/talkowski/rlc47/src/GRCh37.genome \
  #   -i ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz > \
  #   ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
  #   gzip -f ${WRKDIR}/data/CNV/CNV_DENSITY/${group}.${CNV}.coding.CNV_density.bg
  # done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Sanity-check filtering
while read group eti tier descrip include exclude color; do
  for CNV in CNV DEL DUP; do
    echo -e "${group}.${CNV}"
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | wc -l
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l
  done | paste - - - - | awk -v OFS="\t" '{ print $0, $3+$4 }'
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

# #####Gather distributions of CNV per cohort per phenotype
# for CNV in DEL DUP CNV; do
#   echo -e "\n${CNV}\n"
#   for pheno in CTRL DD SCZ CNCR; do
#     while read cohort; do
#       echo -e "${cohort}\t${pheno}\t${CNV}"
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep ${cohort} | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep ${cohort} | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     done < <( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | \
#       fgrep -v "#" | cut -f4 | cut -f1 -d_ | sort | uniq ) | paste - - - - - - -
#     echo -e ""
#   done
#   for pheno in CTRL DD SCZ CNCR; do
#     for dummy in 1; do
#       echo -e "ALL\t${pheno}\t${CNV}"
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     done | paste - - - - - - -
#   done
# done | awk -v OFS="\t" '{ print $1, $2, $3, $4, $6, $8, $5, $7, $9 }'
# #Get median sizes of all CNVs by class & germline/all
# for CNV in DEL DUP CNV; do
#   echo -e "\n${CNV}\n"
#   for dummy in 1; do
#     for pheno in DD SCZ CTRL; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     for pheno in DD SCZ CTRL; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     for pheno in DD SCZ CTRL; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#   done | paste - - -
#   for dummy in 1; do
#     for pheno in DD SCZ CTRL CNCR; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     for pheno in DD SCZ CTRL CNCR; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.coding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     for pheno in DD SCZ CTRL CNCR; do
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.noncoding.bed.gz | fgrep -v "#" | awk '{ print $3-$2 }'
#     done | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#   done | paste - - -
# done

# #####Get reverse CDF of CNV sizes by phenotype
# for group in CTRL DD SCZ DD_SCZ CNCR; do
#   for CNV in DEL DUP; do
#     for dummy in 1; do
#       echo -e "${group}.${CNV}"
#       total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | wc -l )
#       for suffix in 000 0000 00000; do
#         for prefix in $( seq 10 1 90 ); do
#           filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | \
#           awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
#           echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
#         done
#       done | paste -s
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     done | paste -s
#   done
# done > ${WRKDIR}/data/plot_data/CNV_size_bySize.all.txt
# for group in CTRL DD SCZ DD_SCZ CNCR; do
#   for CNV in DEL DUP; do
#     for dummy in 1; do
#       echo -e "${group}.${CNV}"
#       total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | wc -l )
#       for suffix in 000 0000 00000; do
#         for prefix in $( seq 10 1 90 ); do
#           filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | \
#           awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
#           echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
#         done
#       done | paste -s
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.coding.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     done | paste -s
#   done
# done > ${WRKDIR}/data/plot_data/CNV_size_bySize.coding.txt
# for group in CTRL DD SCZ DD_SCZ CNCR; do
#   for CNV in DEL DUP; do
#     for dummy in 1; do
#       echo -e "${group}.${CNV}"
#       total=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | wc -l )
#       for suffix in 000 0000 00000; do
#         for prefix in $( seq 10 1 90 ); do
#           filt=$( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
#           awk -v size="${prefix}${suffix}" '{ if ($3-$2>=size) print $0 }' | wc -l )
#           echo -e "${filt}\t${total}" | awk '{ print $1/$2 }'
#         done
#       done | paste -s
#       zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.noMaxSize.GRCh37.noncoding.bed.gz | fgrep -v "#" | \
#       awk '{ print $3-$2 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#     done | paste -s
#   done
# done > ${WRKDIR}/data/plot_data/CNV_size_bySize.noncoding.txt

# #####Get pct of reference genome covered by at least one CNV
# gsize=$( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | \
#   awk -v OFS="\t" '{ print $1, "1", $2 }' | bedtools subtract -a - \
#   -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed | \
#   awk '{ sum+=($3-$2) }END{ print sum }' )
# #All subjects
# zcat ${WRKDIR}/data/CNV/CNV_MASTER/*CNV.GRCh37.bed.gz | fgrep -v "#" | \
# sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
# -b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
#   '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
# awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
# #Controls only
# zcat ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.CNV.GRCh37.bed.gz | fgrep -v "#" | \
# sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
# -b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
#   '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
# awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
# #Cases only
# zcat ${WRKDIR}/data/CNV/CNV_MASTER/DD_SCZ_CNCR.CNV.GRCh37.bed.gz | fgrep -v "#" | \
# sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
# -b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
#   '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
# awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
# #Cancer only
# zcat ${WRKDIR}/data/CNV/CNV_MASTER/CNCR.CNV.GRCh37.bed.gz | fgrep -v "#" | \
# sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
# -b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
#   '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
# awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'
# #Germline cases only
# zcat ${WRKDIR}/data/CNV/CNV_MASTER/DD.CNV.GRCh37.bed.gz \
# ${WRKDIR}/data/CNV/CNV_MASTER/SCZ.CNV.GRCh37.bed.gz | fgrep -v "#" | \
# sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - | bedtools coverage -a - \
# -b <( grep -ve 'X\|Y\|M' /data/talkowski/rlc47/src/GRCh37.genome | awk -v OFS="\t" \
#   '{ print $1, "1", $2 }' | bedtools subtract -a - -b /data/talkowski/rlc47/src/GRCh37_Nmask.bed ) | \
# awk -v gsize=${gsize} '{ sum+=$5 }END{ print sum/gsize }'

# #####Get files of rCNV sizes for t-test between cases and controls
# mkdir ${WRKDIR}/data/plot_data/CNV_sizes
# for group in CTRL DD SCZ CNCR; do
#   for CNV in DEL DUP; do
#     zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz | \
#     awk '{ print $3-$2 }' > ${WRKDIR}/data/plot_data/CNV_sizes/${group}.${CNV}.sizes.txt
#   done
# done

# #####Test all CNVs for enrichment in genes
# echo -e "all_genes\t${WRKDIR}/data/unfiltered_annotations/gencode.GRCh37.basic_gene_symbols.boundaries.bed.gz" > \
# ${TMPDIR}/elements.txt
# for group in DD_SCZ CNCR; do
#   ${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.CNV.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.CNV.GRCh37.bed.gz \
#   ${TMPDIR}/elements.txt \
#   ${TMPDIR}/${group}_allGenes_enrichment.txt
# done
# #Print results
# for group in DD_SCZ CNCR; do
#   cat ${TMPDIR}/${group}_allGenes_enrichment.txt
# done

#####Run TBRden pileup for all tissue types and all phenotypes (coding & noncoding CNVs)
#Create master mask of N-masked regions and 1Mb flanking telomeres/centromeres
grep -e 'centromere\|telomere' /data/talkowski/rlc47/src/GRCh37_heterochromatin.bed | \
awk -v OFS="\t" '{ print $1, $2-1000000, $3+1000000 }' | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' | \
cat - /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
<( grep -e 'X\|Y\|M' ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed 
#Run TBRden pileups (5kb bins, 5kb step, 1Mb flank dist; run at 0kb, 50kb, and 100kb smoothing)
while read group eti tier descrip include exclude color; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${group} ]; then
    rm -r ${WRKDIR}/analysis/BIN_CNV_pileups/${group}
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${group}
  for smooth in 0 5 10; do
    if ! [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed ]; then
      mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed
    fi
    for CNV in DEL DUP CNV; do
      #Parallelize intersections (LSF)
      bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
      "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 5000 -s 5000 -d 1000000 -r ${smooth} \
      -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.bed \
      -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      /data/talkowski/rlc47/src/GRCh37.genome"
      bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_noncoding \
      "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 5000 -s 5000 -d 1000000 -r ${smooth} \
      -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed \
      -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed  \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
      /data/talkowski/rlc47/src/GRCh37.genome"
      bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_coding \
      "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 5000 -s 5000 -d 1000000 -r ${smooth} \
      -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.coding.bed \
      -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed  \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz \
      /data/talkowski/rlc47/src/GRCh37.genome"
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )
#Rename output directories
while read group eti tier descrip include exclude color; do
  for smooth in 5 10; do
    mv ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}0kb_smoothed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Run TBRden analysis
while read group eti tier descrip include exclude color; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    #Parallelize analyses (LSF)
    bsub -q short -sla miket_sc -J ${group}_${CNV}_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/0kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_all 0.0000003818455 ${color}"
    bsub -q short -sla miket_sc -J ${group}_${CNV}_coding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/0kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.coding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.coding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_coding 0.0000003818455 ${color}"
    bsub -q short -sla miket_sc -J ${group}_${CNV}_noncoding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/0kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_noncoding 0.0000003818455 ${color}"
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" )
#TEST
while read group eti tier descrip include exclude color; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  CNV=DEL
  ${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
  ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/0kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
  ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
  ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
  ${group}_vs_CTRL_${CNV}_noncoding 0.0000003818455 ${color}
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" | fgrep -w SEIZ )

# #Gather example burden files for Jake
# mkdir ${TMPDIR}/example_burden_files
# for smooth in 0 50 100; do
#   mkdir ${TMPDIR}/example_burden_files/${smooth}kb_smoothed
#   for group in CTRL CNCR GERM SEIZ CLNG ASD; do
#     for CNV in CNV DEL DUP; do
#       cp ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.bed.gz \
#       ${TMPDIR}/example_burden_files/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.bed.gz
#       for filt in coding noncoding; do
#         cp ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/0kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz \
#         ${TMPDIR}/example_burden_files/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.bed.${filt}.gz
#       done
#     done
#   done
# done

#####Run 1k CNV shift direct permutation tests for all comparisons
#Note: changed from old 100k matched Fisher permutation
#Note: initial p-value cutoff used: 0.05/130942.8 = 3.818461e-07
#This corresponds to the number of non-overlapping 20kb bins we tested (after blacklisting N-mask, etc)
#20kb chosen since minimum CNV size = 20kb, so max # independent tests = size of genome / 20kb
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
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_all \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
    #Coding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -c -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_coding \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
    #Noncoding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -n -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_noncoding \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
  done
done
# #OLD CODE USED FOR FISHER PERMUTATION
#Split into partitions of 1k permutations each (x100 per comparison)
# for i in $( seq -w 001 100 ); do
#   echo ${i}
#   #Coding + noncoding CNVs
#   if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i} ]; then
#     mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}
#   fi
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 1000000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_all \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
#   #Coding CNVs only
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_coding \
#   -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
#   #Noncoding CNVs only
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_noncoding \
#   -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
# done


#####Collect results from permutation tests
#Collect results across all permutations per group
# for group in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       for i in $( seq -w 001 100 ); do
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_${filt}.permuted.${i}.bed.gz | \
#         fgrep -v "#" | cut -f7 | paste -s
#       done > ${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt
#       Rscript -e "options(scipen=100);\
#       d <- apply(read.table(\"${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt\",header=F),2,sum);\
#       write.table(data.frame(\"perms_less_sig\"=10000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/10001),\
#         \"${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt\",row.names=F,col.names=T,sep=\"\t\",quote=F)"
#       paste <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/001/${group}_vs_CTRL_${CNV}_${filt}.permuted.001.bed.gz | cut -f1-5 ) \
#       ${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt > \
#       ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed
#     done
#   done
# done
#Compose table of all bins for each comparison
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      cat <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -wb -f 1 -r -a - \
      -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
      cut -f1-17,26 ) \
      <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -v -f 1 -r -a - \
      -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
      awk -v OFS="\t" '{ print $0, "NA" }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
      cat <( paste <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
        head -n1 | awk '{ print "#"$0 }' ) <( echo -e "perm_sim_p" ) ) - > \
      ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
    done
  done
done
#Make master table of all bins with summary stats for each bin & each comparison
for dummy in 1; do
  for group in CTRL DD SCZ CNCR; do
    for CNV in DEL DUP; do
      zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
      fgrep -v "#" | cut -f5 > ${TMPDIR}/${group}.${CNV}.all.CNVcounts.tmp
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
        fgrep -v "#" | cut -f17-18 > ${TMPDIR}/${group}.${CNV}.${filt}.pvals.tmp
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
      #Parallelize:
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_filterMasterBurden \
      "${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file_parallelized.sh \
      ${group} ${CNV} ${filt}"
      # list=`mktemp`
      # echo -e "${group}.${CNV}.${filt}.obs_p" > ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R -t 0.05 \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # rm ${list}
      # #Make list of loci (±40kb merge distance & min 20kb size)
      # ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
      # bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct -d 40000 \
      # -i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | \
      # awk '{ if ($3-$2>20000) print $0 }' > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
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
      -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' ) | wc -l
    done
  done
done | awk '{ print "=\""$1"/260\"" }'
#Print overlaps with syndromic CNV intervals for non-coding CNVs
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo -e "\n${group}\n"
  for CNV in ANY_CNV CNV DEL DUP; do
      bedtools intersect -u \
      -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      bedtools intersect -wa -u -b - \
      -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' )
  done | sort -Vk1,1 -k2,2n -k3,3n | uniq
done

#####Print top-5 most significant noncoding loci per comparison
# for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
for group in DD SCZ CNCR; do
  echo -e "\n${group}\n"
  while read chr start end; do
    col=$( zcat ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz | \
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

#####Data for contrasts between deletion and duplication hotspots (Fig 1e)
if [ -e ${WRKDIR}/data/plot_data/del_vs_dup ]; then
  rm -rf ${WRKDIR}/data/plot_data/del_vs_dup
fi
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

#####Count number of autosomal significant loci that don't overlap with consensus syndromic list
for group in DD SCZ DD_SCZ; do
  zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.ANY_FILTER.perm_signif_bins.bed.gz | \
  fgrep -v "#" | cut -f1-3
done | sort -Vk1,1 -k2,2n -k3,3n | uniq | bedtools intersect -u -wa -b - \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | \
cut -f1-3 | bedtools intersect -u -a - \
-b ${WRKDIR}/data/unfiltered_annotations/CollinsFocal_PathogenicCNVs_allSources.boundaries.bed.gz | wc -l

#####Prepare unfiltered annotation files for enrichment tests
#SFARI gene sets
while read genes; do
  echo ${genes}
  fgrep -wf ${SFARI_ANNO}/genelists/${genes}.genes.list \
  ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
  awk '$4 !~ /\-AS/ { print $0 }' | grep -ve '^GL\|^X\|^Y\|^M' | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/${genes}.boundaries.bed
done < <( l ${SFARI_ANNO}/genelists/*genes.list | awk '{ print $9 }' | sed 's/\.genes\.list//g' | \
  xargs -I {} basename {} | fgrep -wv genes_merged_list | fgrep -wv Landrum2014_clinvar_XLinkedDiseaseAssociated )
for group in ASD EPI ID; do
  for level in HC MC LC; do
    echo -e "${group}_${level}"
    #Count
    fgrep -wf /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/genes/${group}_${level}.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed |  \
    awk '$4 !~ /\-AS/ { print $0 }' | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
    ${WRKDIR}/data/unfiltered_annotations/${genes}.boundaries.bed
  done
done
#GTEx top 1k most highly expressed genes by tissue
while read tissue; do
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed.gz | \
  sort -nrk5,5 | grep -ve '^GL\|^X\|^Y\|^M' | head -n1000 | cut -f4 | sort | uniq | \
  fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed |  awk '$4 !~ /\-AS/ { print $0 }' | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_highExpressors_top1k_${tissue}.boundaries.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
#UCNEs
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Dimitrieva2013_UCNE.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/UCNEs.boundaries.bed
#RepeatMasker by class
for class in LINE SINE LTR Satellite Simple_repeat Low_complexity; do
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/src/GRCh37.${class}.RMSK.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/RepMask_${class}.boundaries.bed
done
#HARs
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Doan2016_HAR.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/HAR.boundaries.bed
#SegDups
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/SegDups.boundaries.bed
#ClinGen Pathogenic CNV Regions
zcat ${SFARI_ANNO}/misc/ClinGen_PathogenicCNVs.bed.gz | \
grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/ClinGen_PathoCNV.boundaries.bed
#GWAS catalog variants
zcat ${SFARI_ANNO}/misc/GWAS_Catalog_variants.merged_simple.bed.gz | \
grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/GWAS_Catalog_variants.boundaries.bed
#ENCODE TF ChIP peaks
while read TF; do
  echo ${TF}
  for level in all 0_250 250_500 500_750 750_1000; do
    zcat ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.${level}.bed.gz | \
    grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
    ${WRKDIR}/data/unfiltered_annotations/ENCODE_${TF}_peaks_${level}.boundaries.bed
  done
done < <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | cut -f4 | sort | uniq | \
  cat <( echo -e "all_TF" ) - )
#GTEx eQTLs
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_eQTLs_${tissue}.boundaries.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_eQTLs_AllTissues.boundaries.bed
#COSMIC variants
zcat ${SFARI_ANNO}/misc/COSMIC_variants.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/COSMIC_variants.boundaries.bed
#Affy6 Probes
zcat ${SFARI_ANNO}/misc/Affy6_Probes.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/Affy6_Probes.boundaries.bed
#Illumina 1MDuo Probes
zcat ${SFARI_ANNO}/misc/Illumina_1MDuo_Probes.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/Illumina_1MDuo_Probes.boundaries.bed
#Schmitt TBRs
while read tissue; do
  echo -e "${tissue}"
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/TBRs_${tissue}.boundaries.bed
done < <( cat <( echo "MERGED" ) /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/TBRs_AllTissues.boundaries.bed
#Super enhancers
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/superEnhancers/cleaned/${tissue}_superEnhancer.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/SuperEnhancers_${tissue}.boundaries.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/superEnhancers/*bed | \
  awk '{ print $9 }' | sed 's/superEnhancers\//\t/g' | sed 's/\.bed/\t/g' | cut -f2 )
#Tissue-specific enhancers
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Enhancers_${tissue}.boundaries.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 ) | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Enhancers_all.boundaries.bed
#Brain enhancer RNAs (count)
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Yao2015_BrainEnhancerRNA.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Brain_eRNAs.boundaries.bed
#FMRP targets
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Ascano2012_FMRP.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/FMRP_targets.boundaries.bed
#FANTOM5 NPC, fetal brain, and adult brain enhancers
for tissue in NPC AdultBrain FetalBrain; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/FANTOM52015_${tissue}_TPM5.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/FANTOM5_${tissue}_enhancers.boundaries.bed
done
#Cotney CHD8 peaks
for tissue in Brain hNSC; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Cotney2015_CHD8_${tissue}_peaks.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/CotneyCHD8_${tissue}_peaks.boundaries.bed
done
#Talkowski 36 rCNV loci (count)
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/rCNVs/talkowski_highQual_rCNVs.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Talkowski_rCNVs_HQ.boundaries.bed
#Claire's new rCNV list
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_CER.bed | \
cut -f1,4-5 | sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Redin_PathogenicCNVs_allSources.boundaries.bed
#Final consensus rCNV list
fgrep -v "#" ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $4, "CNV" }' | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed
bedtools intersect -r -f 0.5 -wa -wb \
-a ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed \
-b ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed > \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.selfIntersect.bed
cut -f4 ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed > ${TMPDIR}/CNV_interval_IDs.list
source /apps/lab/miket/anaconda/4.0.5/envs/collins_py3/bin/activate collins_py3
/data/talkowski/rlc47/code/svcf/scripts/bedcluster -p PathogenicCNVs_allSources_nonredundant_hg19_RLC -m \
${TMPDIR}/CNV_interval_IDs.list \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.selfIntersect.bed \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed
while read CNVID; do
  awk -v OFS="\t" -v CNVID=${CNVID} '{ if ($7==CNVID) print $1, $2, $3, $5 }' \
  ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed | \
  bedtools merge -c 4 -o distinct -i - 
done < <( fgrep -v "#" ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed | \
  cut -f7 | sort -Vk1,1 | uniq ) | sort -Vk1,1 -k2,2n -k3,3n > \
${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_RLC.bed
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_RLC.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/CollinsFocal_PathogenicCNVs_allSources.boundaries.bed
#rCNV list by source
while read study; do
  echo ${study}
  awk -v study=${study} '{ if ($4==study) print $0 }' \
  ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
  sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/unfiltered_annotations/${study}_PathoCNVs.boundaries.bed
done < <( fgrep -v "#" ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
  cut -f4 | sort | uniq )
#Gzip all annotations
gzip -f ${WRKDIR}/data/unfiltered_annotations/*.boundaries.bed

#####Launch unfiltered hotspot functional enrichment tests (both binomial and t-test)
if ! [ -e ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered ]; then
  mkdir ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered
fi
while read anno; do
  if [ -e ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/ ]; then
    rm -rf ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/
  fi
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/
  for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
    for CNV in ANY_CNV CNV DEL DUP; do
      for filt in ANY_FILTER all coding noncoding; do 
        #Binomial test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t binomial \
        -p ${group}_${CNV}_${filt}.${anno}.unfiltered_binomial \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/"
        #t test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t t \
        -p ${group}_${CNV}_${filt}.${anno}.unfiltered_t \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/"
      done
    done
  done > ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh
  chmod a+x ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh
  bsub -q short -sla miket_sc -J hotspot_unfiltered_enrichment_${anno} -u nobody \
  "sh ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh"
# done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
#         awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )
#optional: exclude ENCODE, GTEx, Enhancer, Super Enhancer, and TAD datasets
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve 'GTEx\|ENCODE\|Enhancer\|TBRs_' )

#####Launch raw CNV set enrichment tests
paste <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve '250_500\|500_750' ) \
<( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | grep -ve '250_500\|500_750' ) > \
${TMPDIR}/annotations_to_test.list
#Run first pass with adjustment to size
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/CNV_set_enrichments/${group} ]; then
    rm -rf ${WRKDIR}/analysis/CNV_set_enrichments/${group}
  fi
  mkdir ${WRKDIR}/analysis/CNV_set_enrichments/${group}
  for CNV in CNV DEL DUP; do
    #All
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt"
    for filt in coding noncoding; do 
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt"
    done
  done
done
#Run second pass with adjustment to all genes, not to size
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    #All
    adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt.gz | \
      fgrep -w gencode.GRCh37.basic_gene_symbols | cut -f6 )
    echo ${adjust}
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment_geneNorm -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.geneNorm.setEnrichments.txt"
    for filt in coding noncoding; do 
      adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt.gz | \
        fgrep -w gencode.GRCh37.basic_gene_symbols | cut -f6 )
      echo ${adjust}
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment_geneNorm -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.geneNorm.setEnrichments.txt"
    done
  done
done
#Run third pass with adjustment based on median enrichment for all annotations tested
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    #All
    adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt.gz | \
      cut -f6 | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
    echo ${adjust}
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment_medEnrichNorm -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.medEnrichNorm.setEnrichments.txt"
    for filt in coding noncoding; do 
      adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt.gz | \
        cut -f6 | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
      echo ${adjust}
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment_medEnrichNorm -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.medEnrichNorm.setEnrichments.txt"
    done
  done
done


#####Filter all annotations by excluding elements within 10kb of a protein-coding exon
while read anno; do
  echo ${anno}
  bedtools intersect -v -wa \
  -a ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed > \
  ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed
  if [ -s ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed ]; then
    gzip -f ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed 
  else
    echo "Removing zero-sized annotation file for ${anno}"
    rm ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed
  fi
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | xargs -I {} basename {} | sed 's/\.boundaries\.bed\.gz//g' )

#####Launch filtered hotspot functional enrichment tests (both binomial and t-test) for noncoding sites only
if ! [ -e ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered ]; then
  mkdir ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered
fi
while read anno; do
  if [ -e ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/ ]; then
    rm -rf ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/
  fi
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/
  for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
    for CNV in ANY_CNV CNV DEL DUP; do
      for filt in noncoding; do 
        #Binomial test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t binomial \
        -p ${group}_${CNV}_${filt}.${anno}.filtered_binomial \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/"
        #t test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t t \
        -p ${group}_${CNV}_${filt}.${anno}.filtered_t \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/"
      done
    done
  done > ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh
  chmod a+x ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh
  bsub -q short -sla miket_sc -J hotspot_filtered_enrichment_${anno} -u nobody \
  "sh ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh"
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve '250_500\|500_750' )

#####Collect results from unfiltered enrichment tests
if ! [ -e ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS ]; then
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS
fi
#Print headers to output files
for dummy in 1; do
  for prefix in binomial_OR binomial_CI_min binomial_CI_max binomial_p \
  t_fold t_CI_min t_CI_max t_p; do
    for dummy in 2; do
      echo "annotation"
      for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
        for CNV in ANY_CNV CNV DEL DUP; do
          for filt in ANY_FILTER all coding noncoding; do
            echo -e "${group}_${CNV}_${filt}"
          done
        done
      done | paste -s
    done | paste - - > ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.${prefix}.txt
  done
done
while read anno; do
  echo ${anno}
  #Binomial odds ratio
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_OR.txt
  #Binomial odds ratio CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_CI_min.txt
  #Binomial odds ratio CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_CI_max.txt
  #Binomial p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_p.txt
  #t-test fold-enrichment
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_fold.txt
  #t-test CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_CI_min.txt
  #t-test CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_CI_max.txt
  #t-test p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_p.txt
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )

#####Collect results from filtered enrichment tests
if ! [ -e ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS ]; then
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS
fi
#Print headers to output files
for dummy in 1; do
  for prefix in binomial_OR binomial_CI_min binomial_CI_max binomial_p \
  t_fold t_CI_min t_CI_max t_p; do
    for dummy in 2; do
      echo "annotation"
      for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
        for CNV in ANY_CNV CNV DEL DUP; do
          for filt in ANY_FILTER all coding noncoding; do
            echo -e "${group}_${CNV}_${filt}"
          done
        done
      done | paste -s
    done | paste - - > ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.${prefix}.txt
  done
done
while read anno; do
  echo ${anno}
  #Binomial odds ratio
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_OR.txt
  #Binomial odds ratio CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_CI_min.txt
  #Binomial odds ratio CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_CI_max.txt
  #Binomial p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_p.txt
  #t-test fold-enrichment
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_fold.txt
  #t-test CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_CI_min.txt
  #t-test CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_CI_max.txt
  #t-test p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_p.txt
done < <( l ${WRKDIR}/data/filtered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )

#####Gather data for published CNV locus enrichments (Fig1f)
#Get ORs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f4
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f4
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt
#Get lower CIs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.lowerCI.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f5
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f5
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.lowerCI.txt
#Get upper CIs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.upperCI.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f6
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f6
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.upperCI.txt

#####Prepare annotation files for functional enrichment analyses
#Make master bin file
zcat ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz | \
cut -f1-4 | sed '1d' > ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
#Subset master bin file to autosomes only
zcat ${WRKDIR}/data/annotations/GRCh37.master_bins.bed.gz | fgrep -v "X" | fgrep -v "Y" | sed '1d' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
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
#PhyloP 100-way conservation
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
/scratch/miket/rlc47temp/DOSAGE/hg19.100way.phyloP100way.bw \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
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

#####Get noncoding TBR fold-enrichments for noncoding CNVs
for group in ANY_DISEASE DD_SCZ CNCR; do
  echo ${group}
  for CNV in CNV DEL DUP; do
    tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/TBRs_MERGED/${group}_${CNV}_noncoding.TBRs_MERGED.filtered_binomial.TBRden_binomial_annotation_test.results.txt | \
    awk '{ print $2/$3 }'
    tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/TBRs_MERGED/${group}_${CNV}_noncoding.TBRs_MERGED.filtered_binomial.TBRden_binomial_annotation_test.results.txt | \
    cut -f5,6
  done | paste -s 
done | paste - - > ${WRKDIR}/data/plot_data/Noncoding_rCNV_hotspot.noncoding_TBR_enrichments.txt

#####Run protein-coding exon pileups for noncoding rCNVs in flanks
#Prep list of 40kb flanks
distance=40000
awk -v OFS="\t" -v d=${distance} '{ print $1, $2-d, $3+d, $4 }' \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed | \
awk -v OFS="\t" '{ if ($2<0) $2=0; if ($3<1) $3=1; print }' > \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.40kb_windows.bed
#Launch tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/ ]; then
      rm -rf ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/
    fi
    mkdir ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_exon_flanking_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -n -z -t upper \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    -p ${group}_${CNV}_noncoding_exonFlanks \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.40kb_windows.bed \
    ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/"
  done
done

#####Run TBRden TBR CNV pileups
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/ ]; then
      rm -rf ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/
    fi
    mkdir ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBR_pileup_all \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_pileup.sh -z \
    -o ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_all_TBR_pileups.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBR_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_pileup.sh -z \
    -o ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_noncoding_TBR_pileups.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz"
  done
done

#####Run TBRden TBR CNV burden tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ ]; then
      rm -rf ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/
    fi
    mkdir ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/
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
    #Parallelize (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_all \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/TBR_CNV_pileups/CTRL_${CNV}_TBR_pileups/CTRL_${CNV}_all_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_all_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ \
    ${group}_${CNV}_all 0.000009072764 ${color}"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/TBR_CNV_pileups/CTRL_${CNV}_TBR_pileups/CTRL_${CNV}_noncoding_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_noncoding_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ \
    ${group}_${CNV}_noncoding 0.000009072764 ${color}"
  done
done

#####Run direct permutation tests of TBRs 
#Launch tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_TBR_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -n -z -t upper \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    -p ${group}_${CNV}_TBR_burdens_noncoding \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_TBR_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -z -t upper \
    -p ${group}_${CNV}_TBR_burdens_all \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/"
  done
done

#####Aggregate results of TBR noncoding rCNV burden tests
#Print header
for dummy in 1; do
  echo -e "#chr\tstart\tend\tTBR_ID\ttissues"
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in DEL DUP CNV; do
      echo -e "${group}_${CNV}_p_obs\n${group}_${CNV}_p_perm"
    done
  done | paste -s
done | paste -s > ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
#Gather p-values
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    zcat ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/${group}_${CNV}_noncoding.TBRden_results.bed.gz | \
    sed '1d' | awk '{ print $NF }' | paste -s
    zcat ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/${group}_${CNV}_TBR_burdens_noncoding.TBRden_direct_test_results.bed.gz | \
    sed '1d' | awk '{ print $(NF-1) }' | paste -s
  done
done > ${TMPDIR}/matrix_pretranspose.txt
#Transpose matrix
Rscript -e "write.table(t(read.table(\"${TMPDIR}/matrix_pretranspose.txt\",header=F)),\"${TMPDIR}/matrix_posttranspose.txt\",col.names=F,row.names=F,sep=\"\\t\")"
paste <( zcat ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz | fgrep -v "#" | cut -f1-5 ) \
${TMPDIR}/matrix_posttranspose.txt >> ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
gzip -f ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
#Copy to plot directory
cp ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt.gz ${WRKDIR}/data/plot_data/



