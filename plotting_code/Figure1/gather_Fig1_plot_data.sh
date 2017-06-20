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
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

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

####Get number of phenotypic group assignments per patient
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | \
sed 's/,/\t/g' | awk '{ print NF }' | sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | \
sort -nk1,1 > ${WRKDIR}/data/plot_data/figure1/phenotypes_per_sample_hist.txt

####Get count of patients with a NEURO and/or a SOMA HPO assignment
#Both
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | fgrep -v TCGA | \
fgrep -wf <( fgrep -w "NEURO" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -wf \
<( fgrep -w "SOMA" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -v "0002664" | wc -l
#NEURO, no SOMA
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | fgrep -v TCGA | \
fgrep -wf <( fgrep -w "NEURO" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -wvf \
<( fgrep -w "SOMA" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -v "0002664" | wc -l
#SOMA, no NEURO
cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | fgrep -v TCGA | \
fgrep -wvf <( fgrep -w "NEURO" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -wf \
<( fgrep -w "SOMA" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' ) | fgrep -v "0002664" | wc -l

#####Get matrix of count of patients overlapping per phenotype group
while read dA; do
  while read groupA etiA tierA descriptionA includeA excludeA colorA nA; do
    echo ${groupA}
    while read dB; do
      while read groupB etiB tierB descriptionB includeB excludeB colorB nB; do
        if [ ${groupA} == ${groupB} ]; then
          echo ${nA}
        else
          if [ ${etiA} != ${etiB} ]; then
            echo 0
          else
            cut -f2 ${WRKDIR}/data/HPO_map/master_patient_IDs_and_phenos.wHPO.list | \
            fgrep -wf <( echo ${includeA} | sed 's/\;/\n/g' ) | \
            fgrep -wvf <( echo ${excludeA} | sed 's/\;/\n/g' ) | \
            fgrep -wf <( echo ${includeB} | sed 's/\;/\n/g' ) | \
            fgrep -wvf <( echo ${excludeB} | sed 's/\;/\n/g' ) | wc -l
          fi
        fi
      done < <( fgrep -w ${dB} ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )
    done < <( echo -e "Neurological_defect\nNeurodevelopmental_disorder\n\
Developmental_delay\nNeuropsychiatric_disorder\nSchizophrenia\nAutism_spectrum_disorder\n\
Seizures\nHypotonia\nBehavioral_abnormality_other\nIntellectual_disability\n\
Nonneurological_defect\nHead_or_neck_defect\nGrowth_defect\nCardiac_defect\n\
Skeletal_defect\nDigestive_respiratory_or_urinary_defect\nMuscular_defect\n\
Eye_or_ear_defect\nIntegument_defect\nEndocrine_metabolic_or_immune_defect" ) 
  done < <( fgrep -w "${dA}" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | paste -s
done < <( echo -e "Neurological_defect\nNeurodevelopmental_disorder\n\
Developmental_delay\nNeuropsychiatric_disorder\nSchizophrenia\nAutism_spectrum_disorder\n\
Seizures\nHypotonia\nBehavioral_abnormality_other\nIntellectual_disability\n\
Nonneurological_defect\nHead_or_neck_defect\nGrowth_defect\nCardiac_defect\n\
Skeletal_defect\nDigestive_respiratory_or_urinary_defect\nMuscular_defect\n\
Eye_or_ear_defect\nIntegument_defect\nEndocrine_metabolic_or_immune_defect" ) > \
${WRKDIR}/data/plot_data/figure1/germline_case_overlap.matrix.txt

#####Get sizes of all CNVs and all filters per CNV class and tier 2 phenotype group
for group in CTRL NEURO NDD PSYCH SOMA CNCR; do
  for filt in all coding haplosufficient noncoding intergenic; do
    for VF in E2 E3 E4 N1; do
      for CNV in CNV DEL DUP; do
        zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.E2.GRCh37.${filt}.bed.gz | \
        fgrep -v "#" | awk '{ print $3-$2 }' > \
        ${WRKDIR}/data/plot_data/figure1/${CNV}_size.${group}.noMaxSize.${freq}.${filt}.txt
      done
    done
  done
done
#Binned by VF - PER VARIANT
for group in CTRL NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    for max in E2 E3 E4; do
      for min in E3 E4 N1; do
        zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.${max}.GRCh37.all.bed.gz | \
        fgrep -v "#" | fgrep -wvf \
        <( zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.${min}.GRCh37.all.bed.gz | \
        fgrep -v "#" | cut -f4 ) | cut -f4 | fgrep -wf - \
        ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
        awk -v OFS="\t" '{ print $3-$2, $7 }' | sort -nk1,1 -k2,2 | uniq > \
        ${WRKDIR}/data/plot_data/figure1/${CNV}_size.${group}.noMaxSize.${max}_${min}.all.byVariant.txt
      done
    done
    #N1
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.noMaxSize.${max}.GRCh37.all.bed.gz | \
    fgrep -v "#" | cut -f4 | fgrep -wf - \
    ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
    awk -v OFS="\t" '{ print $3-$2, $7 }' | sort -nk1,1 -k2,2 | uniq > \
    ${WRKDIR}/data/plot_data/figure1/${CNV}_size.${group}.noMaxSize.N1_N1.all.byVariant.txt
  done
done
#CNCR
for CNV in CNV DEL DUP; do
  for VF in E2_E3 E3_E4 E4_N1 N1_N1; do
    cp ${WRKDIR}/data/plot_data/figure1/${CNV}_size.CNCR.noMaxSize.E2.all.txt \
    ${WRKDIR}/data/plot_data/figure1/${CNV}_size.CNCR.noMaxSize.${VF}.all.byVariant.txt
  done
done


#####Get sizes and VFs of all CNVs per tier 2 for scatterplots -- nixed
#Germline CNVs
for group in CTRL NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    zcat ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.E2.GRCh37.all.bed.gz | \
    fgrep -v "#" | cut -f4 | fgrep -wf - \
    ${WRKDIR}/data/CNV/CNV_RAW/merged_CNV/germline_${CNV}.merged.bed | \
    cut -f2,3,7,9 | uniq | awk -v OFS="\t" '{ print $2-$1, $4/102257 }' > \
    ${WRKDIR}/data/plot_data/figure1/${CNV}_size_and_VF.${group}.txt
  done
done











