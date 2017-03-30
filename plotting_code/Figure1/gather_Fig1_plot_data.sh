#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure 1

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

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/figure1 ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure1
fi
mkdir ${WRKDIR}/data/plot_data/figure1

#####Get counts of patients per analysis group
for dummy in 1; do
  #CTRL
  echo -e "CTRL\t#A5A6A7\tblack\t38628\tHealthy_controls"
  echo -e "SKIP\tNA\tNA\t0\tSKIP"
  #Germline all
  while read group eti tier description include exclude color n; do
    if [ ${group} == "GERM" ]; then
      echo -e "${group}\t#7B2AB3\tblack"
    else
      echo -e "${group}\t#CAAAE1\tblack"
    fi
    if [ ${eti} == "CNCR" ]; then
      echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
      <( fgrep TCGA ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | cut -f2 ) | \
      fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
    else
      if [ ${group} = "CTRL" ]; then
        echo  -e "38628"
      else
        echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
        <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
        fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
      fi
    fi
    echo -e "${description}"
  done < <( grep -e 'Germline_disease\|Unknown' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  paste - - - | sort -nrk4,4
  echo -e "SKIP\tNA\tNA\t0\tSKIP"
  #Neurological
  while read group eti tier description include exclude color n; do
    if [ ${group} == "NEURO" ]; then
      echo -e "${group}\t#00BFF4\tblack"
    else
      echo -e "${group}\t#99E5FB\tblack"
    fi
    if [ ${eti} == "CNCR" ]; then
      echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
      <( fgrep TCGA ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | cut -f2 ) | \
      fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
    else
      if [ ${group} = "CTRL" ]; then
        echo  -e "38628"
      else
        echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
        <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
        fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
      fi
    fi
    echo -e "${description}"
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
    grep -e 'Neurological_defect\|Neurodevelopmental_disorder\|Autism_spectrum_disorder\|Developmental_delay\|Intellectual_disability\|Neuropsychiatric_disorder\|Schizophrenia\|Behavioral_abnormality_other\|Seizures\|Hypotonia' ) | \
  paste - - - | sort -nrk4,4
  echo -e "SKIP\tNA\tNA\t0\tSKIP"
  #Non-neurological
  while read group eti tier description include exclude color n; do
    if [ ${group} == "SOMA" ]; then
      echo -e "${group}\t#EC008D\tblack"
    else
      echo -e "${group}\t#F799D1\tblack"
    fi
    if [ ${eti} == "CNCR" ]; then
      echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
      <( fgrep TCGA ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | cut -f2 ) | \
      fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
    else
      if [ ${group} = "CTRL" ]; then
        echo  -e "38628"
      else
        echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
        <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
        fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
      fi
    fi
    echo -e "${description}"
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
    grep -e 'Nonneurological_defect\|Digestive_respiratory_or_urinary_defect\|Growth_defect\|Integument_defect\|Skeletal_defect\|Head_or_neck_defect\|Muscular_defect\|Cardiac_defect\|Eye_or_ear_defect\|Endocrine_metabolic_or_immune_defect' ) | \
  paste - - - | sort -nrk4,4
  echo -e "SKIP\tNA\tNA\t0\tSKIP"
  #Cancer
  while read group eti tier description include exclude color n; do
    if [ ${group} == "CNCR" ]; then
      echo -e "${group}\t#FFCB00\tblack"
    else
      echo -e "${group}\t#FFEA99\tblack"
    fi
    if [ ${eti} == "CNCR" ]; then
      echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
      <( fgrep TCGA ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | cut -f2 ) | \
      fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
    else
      if [ ${group} = "CTRL" ]; then
        echo  -e "38628"
      else
        echo -e "${include}" | sed 's/\;/\n/g' | fgrep -wf - \
        <( cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list ) | \
        fgrep -wvf <( echo -e "${exclude}" | sed 's/\;/\n/g' ) | wc -l
      fi
    fi
    echo -e "${description}"
  done < <( fgrep CNCR ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | paste - - - | sort -nrk4,4
done > ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt

#####Get sizes of all E2 CNVs per filter tier 2 phenotype group
for group in CTRL NEURO SOMA CNCR; do
  for filt in coding haplosufficient noncoding intergenic; do
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.CNV.noMaxSize.E2.GRCh37.${filt}.bed.gz | \
    fgrep -v "#" | awk '{ print $3-$2 }' > \
    ${WRKDIR}/data/plot_data/figure1/CNV_size.${group}.noMaxSize.E2.${filt}.txt
  done
done

#####Get sizes of all CNVs per VF filter per tier 2 phenotype group
for group in CTRL NEURO SOMA CNCR; do
  for freq in E2 E3 E4 N1; do
    for CNV in DEL DUP; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.${freq}.GRCh37.all.bed.gz | \
      fgrep -v "#" | awk '{ print $3-$2 }' > \
      ${WRKDIR}/data/plot_data/figure1/${CNV}_size.${group}.${freq}.all.txt
    done
  done
done












