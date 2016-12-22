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
TMPDIR=${WRKDIR}/tmp
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
mkdir ${WRKDIR}/data/CNV/CNV_MASTER
mkdir ${WRKDIR}/data/plot_data
mkdir ${WRKDIR}/data/annotations
mkdir ${WRKDIR}/data/annotations/exonExclusion
mkdir ${WRKDIR}/analysis
mkdir ${WRKDIR}/analysis/BIN_CNV_pileups
mkdir ${WRKDIR}/analysis/BIN_CNV_burdens
mkdir ${WRKDIR}/analysis/BIN_CNV_permutation
mkdir ${WRKDIR}/analysis/Final_Loci
mkdir ${WRKDIR}/analysis/Functional_Enrichments
mkdir ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion

#####Copy CNVs from TAD intolerance Project
#Note: see old code for parameters & syntax used during callset curation
cp /data/talkowski/rlc47/TAD_intolerance/data/CNV/CNV_MASTER/* \
${WRKDIR}/data/CNV/CNV_MASTER/

#####Run TBRden pileup for all tissue types and all phenotypes (coding & noncoding CNVs)
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
  done
done

#####Run TBRden analysis
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    #Parallelize analyses (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_all"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_noncoding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_noncoding"
  done
done

#####Run 10k CNV shift permutation tests for all comparisons
#Note: initial p-value cutoff used: 0.05/26802 = 1.865532e-06
#This corresponds to the number of non-overlapping autosomal 100kb bins we tested (after blacklisting N-mask, etc)
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    #Coding + noncoding CNVs
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
    bsub -q normal -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.10k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 300000 -N 10000 \
    -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.permuted.bed \
    -p ${group}_vs_CTRL_${CNV}_all \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
    #Noncoding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
    bsub -q normal -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.10k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 300000 -N 10000 \
    -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.permuted.bed \
    -p ${group}_vs_CTRL_${CNV}_noncoding \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.noncoding.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
  done
done


