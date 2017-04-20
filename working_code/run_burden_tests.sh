#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master annotation curation code
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
mkdir ${WRKDIR}/analysis/annoSet_burden
while read pheno; do
  if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}
    for filt in all coding haplosufficient noncoding intergenic; do
      if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt} ]; then
        rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}
      fi
      mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}
      for VF in E2 E3 E4 N1; do
        if [ -e ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF} ]; then
          rm -rf ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}
        fi
        mkdir ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Run annotation set burden testing
#iterate & prep directories
while read pheno; do
  for CNV in CNV DEL DUP; do
    #CODE:
    # for filt in all coding haplosufficient noncoding intergenic; do
    #   for VF in E2 E3 E4 N1; do
    #     #Submit batch job
    #     ${WRKDIR}/bin/rCNVmap/bin/annoSet_burdenTest_batch.sh -N 1000 \
    #       -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    #       -p ${pheno}_${CNV}_${filt}_${VF} \
    #       -o ${WRKDIR}/analysis/annoSet_burden/${pheno}/${CNV}/${filt}/${VF}/ \
    #       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
    #       ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
    #       ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list \
    #       ${WRKDIR}/data/misc/GRCh37_autosomes.genome
    #   done
    # done
    #PARALLELIZE:
    bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_annoSet_burdens \
    ""
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )




