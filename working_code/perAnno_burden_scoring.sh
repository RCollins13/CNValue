#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run all rCNV burden scoring per annotation

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directories if exists
if [ -e ${WRKDIR}/data/perAnno_burden ]; then
  rm -rf ${WRKDIR}/data/perAnno_burden
fi
mkdir ${WRKDIR}/data/perAnno_burden
if [ -e ${WRKDIR}/analysis/perAnno_burden ]; then
  rm -rf ${WRKDIR}/analysis/perAnno_burden
fi
mkdir ${WRKDIR}/analysis/perAnno_burden

#####Create subdirectories
while read pheno; do
  if [ -e ${WRKDIR}/data/perAnno_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perAnno_burden/${pheno}
  if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Reorder noncoding annotations based on number of elements in group
while read class path; do
  nElements=$( cat ${path} | wc -l )
  echo -e "${class}\t${path}\t${nElements}"
done < ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list | \
sort -nk3,3 | cut -f1-2 > \
${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list

#####Submit burden data collection for all phenotypes
while read pheno; do
  for CNV in CNV DEL DUP; do
    # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV} ]; then
    #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}
    # fi
    # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}
    for VF in E4; do
      # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF} ]; then
      #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}
      # fi
      # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}
      for filt in haplosufficient noncoding; do
        # if [ -e ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt} ]; then
        #   rm -rf ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        # fi
        # mkdir ${WRKDIR}/data/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perAnno_burden_dataCollection_${filt} \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/gather_annoScore_data_batchMode.sh -F \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Submit annoScore model for all phenotypes
while read pheno; do
  #Get number of subjects in group
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
           ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}
    for VF in E4; do
      if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF} ]; then
        rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}
      fi
      mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}
      for filt in haplosufficient noncoding; do
        if [ -e ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt} ]; then
          rm -rf ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        fi
        mkdir ${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}
        #Launch model in batch mode across all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_annoScoreModel \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/run_annoScore_model_batchMode.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect significant elements per track per phenotype
while read pheno; do
  #Make output directories (if necessary)
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}
  fi
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL
  fi
  for CNV in CNV DEL DUP; do
    for VF in E4; do
      for filt in haplosufficient noncoding; do
        #Launch collection script for all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_annoScoreModel \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_perClass_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.prioritized_for_annoScore_modeling.list"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Collect significant elements across all tracks per phenotype
#Make lists of tissue-defined elements
while read tissue; do
  fgrep ${tissue} ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list > \
  ${WRKDIR}/lists/all_${tissue}_genome_annotations.list
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
          sort | uniq )
#Make list of tissue-agnostic annotations
cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
sort | uniq | fgrep -vf - \
${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list > \
${WRKDIR}/lists/tissue_agnostic_genome_annotations.list
#Iterate over phenotypes
while read pheno; do
  #Make output directories (if necessary)
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged
  fi
  for CNV in CNV DEL DUP; do
    for VF in E4; do
      for filt in haplosufficient noncoding; do
        #Launch merge & filtering script across all annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_all \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list \
        bonf \
        ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.all_classes.bed"
        #Launch merge & filtering script for tissue-agnostic annotations
        bsub -q normal -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_tissueAgnostic \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
        ${pheno} ${CNV} ${VF} ${filt} \
        ${WRKDIR}/lists/tissue_agnostic_genome_annotations.list \
        bonf \
        ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.tissue_agnostic.bed"
        #Iterate over tissues and launch merge & filtering script for tissue-dependent annotations
        while read tissue; do
          bsub -q short -u nobody -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_mergeSignificantElements_${tissue} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_significant_elements_merged_perPheno_annoScore.sh \
          ${pheno} ${CNV} ${VF} ${filt} \
          ${WRKDIR}/lists/all_${tissue}_genome_annotations.list \
          bonf \
          ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.${tissue}.bed"
        done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq )
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Get count of significant loci by phenotype
#Merged across all annotations
VF=E4
while read pheno; do
  for dummy in 1; do
    echo ${pheno}
    for CNV in DEL DUP; do
      for filt in haplosufficient noncoding; do
        fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.all_classes.bed | \
        awk '{ if ($3-$2>=5000 && $3-$2<500000) print $0 }' | wc -l
      done
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Merge significant loci across phenotypes & between DEL/DUP
#Create working directory
if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/ ]; then
  mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/
fi
#Bedcluster across phenotypes with 50% reciprocal overlap required
#Min size = 5kb; applied *BEFORE* merging
#Max size = 500kb; applied *BEFORE* merging
minSize=5000
maxSize=500000
VF=E4
while read annoSet; do
  for CNV in DEL DUP; do
    for filt in haplosufficient noncoding; do
      echo -e "${annoSet}_${CNV}_${filt}"
      #Create master list of all significant elements
      while read pheno; do
        fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.${CNV}.${VF}.${filt}.bonf_sig_elements_merged.all_classes.bed | \
        awk -v OFS="\t" -v pheno=${pheno} -v CNV=${CNV} -v filt=${filt} -v minSize=${minSize} -v maxSize=${maxSize} \
        '{ if ($3-$2>=minSize && $3-$2<=maxSize) print $1, $2, $3, pheno"_"CNV"_"filt"_"NR, pheno"_"CNV"_"filt"_"NR, CNV }' 
      done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
                cut -f1 | fgrep -v CTRL ) | sort -Vk1,1 -k2,2n -k3,3n > \
      ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.bed
      #Run bedtools intersect (50% recip)
      bedtools intersect -r -f 0.5 -wa -wb \
      -a ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.bed \
      -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.bed > \
      ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.all_vs_all.bed
      #Run bedcluster
      cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.bed | \
      sort | uniq > ${TMPDIR}/${annoSet}_${CNV}_${VF}_${filt}.input_element_IDs.tmp
      /data/talkowski/rlc47/code/svcf/scripts/bedcluster -p ${annoSet}_${CNV}_${VF}_${filt}_mergedSignificantLoci -m \
      ${TMPDIR}/${annoSet}_${CNV}_${VF}_${filt}.input_element_IDs.tmp \
      ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.pre_merge.all_vs_all.bed \
      ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.merged.bed
    done
  done
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - )
#Bedcluster across DEL/DUP with 50% reciprocal overlap required
while read annoSet; do
  for filt in haplosufficient noncoding; do
    echo -e "${annoSet}_${filt}"
    #Create master list of all significant elements
    for CNV in DEL DUP; do
      fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_${CNV}_${VF}_${filt}.signif_loci.merged.bed | \
      awk -v OFS="\t" '{ print $1, $2, $3, $7, $7, "element" }' | uniq
    done | sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.bed
    #Run bedtools intersect (50% recip)
    bedtools intersect -r -f 0.5 -wa -wb \
    -a ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.bed \
    -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.bed > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.all_vs_all.bed
    #Run bedcluster
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.all_vs_all.bed | \
    sort | uniq > ${TMPDIR}/${annoSet}_DELDUP_${VF}_${filt}.input_element_IDs.tmp
    /data/talkowski/rlc47/code/svcf/scripts/bedcluster -p ${annoSet}_DELDUP_${VF}_${filt}_mergedSignificantLoci -m \
    ${TMPDIR}/${annoSet}_DELDUP_${VF}_${filt}.input_element_IDs.tmp \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.pre_merge.all_vs_all.bed \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DELDUP_${VF}_${filt}.signif_loci.merged.bed
  done
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - )
#Bedcluster for haplosufficient DEL and noncoding DUP (analysis for paper)
while read annoSet; do
  echo -e "${annoSet}"
  #Create master list of all significant elements
  cat <( fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DEL_${VF}_haplosufficient.signif_loci.merged.bed | \
          awk -v OFS="\t" '{ print $1, $2, $3, $7, $7, "element" }' | uniq ) \
      <( fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_DUP_${VF}_noncoding.signif_loci.merged.bed | \
          awk -v OFS="\t" '{ print $1, $2, $3, $7, $7, "element" }' | uniq ) | \
      sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.bed
  #Run bedtools intersect (50% recip)
  bedtools intersect -r -f 0.5 -wa -wb \
  -a ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.bed \
  -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.bed > \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.all_vs_all.bed
  #Run bedcluster
  cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.all_vs_all.bed | \
  sort | uniq > ${TMPDIR}/${annoSet}_haplosuffDELnoncodingDUP_${VF}.input_element_IDs.tmp
  /data/talkowski/rlc47/code/svcf/scripts/bedcluster -p ${annoSet}_haplosuffDELnoncodingDUP_${VF}_mergedSignificantLoci -m \
  ${TMPDIR}/${annoSet}_haplosuffDELnoncodingDUP_${VF}.input_element_IDs.tmp \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.pre_merge.all_vs_all.bed \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.bed
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - )
#If two elements overlap, keep the smaller of the two
VF=E4
while read annoSet; do
  echo -e "${annoSet}"
  #First, run overlap
  bedtools intersect -wa -wb \
  -a <( cut -f1-3,7 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.bed ) \
  -b <( cut -f1-3,7 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.bed ) > \
  ${TMPDIR}/all_overlaps.bed
  #Second, get list of elements that don't overlap a different element
  cut -f4 ${TMPDIR}/all_overlaps.bed | sort -Vk1,1 | uniq -c | \
  awk -v OFS="\t" '{ if ($1==1) print $2 }' > \
  ${TMPDIR}/elements_to_keep.IDs.list
  #Third, iterate over all remaining elements and only keep the smallest overlapping element
  # among all possible overlapping elements per element
  while read ID; do
    awk -v ID=${ID} -v OFS="\t" '{ if ($4==ID || $8==ID) print $4, $3-$2"\n"$8, $7-$6 }' \
    ${TMPDIR}/all_overlaps.bed | sort -nk2,2 | head -n1 | cut -f1
  done < <( cut -f4 ${TMPDIR}/all_overlaps.bed | sort -Vk1,1 | uniq | \
            fgrep -wvf ${TMPDIR}/elements_to_keep.IDs.list ) | \
  sort -Vk1,1 | uniq >> ${TMPDIR}/elements_to_keep.IDs.list
  #Finally, get those IDs from the original bed file and write to new file
  fgrep -wf ${TMPDIR}/elements_to_keep.IDs.list \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.bed | \
  cut -f1-3,7 | sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | uniq > \
  ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed
done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - )

#Split final clustered loci by phenotype
#Create output directory
if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci ]; then
  mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci
fi
while read pheno; do
  echo ${pheno}
  if [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}
  while read annoSet; do
    #Haplosuff DEL
    bedtools intersect -f 0.5 -a ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed \
    -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.DEL.${VF}.haplosufficient.bonf_sig_elements_merged.all_classes.bed | \
    sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | uniq > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.${annoSet}.bed
    #Noncoding DUP
    bedtools intersect -f 0.5 -a ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed \
    -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/merged/${pheno}.DUP.${VF}.noncoding.bonf_sig_elements_merged.all_classes.bed | \
    sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | uniq > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.${annoSet}.bed
  done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
                  sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - )
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Get counts of significant loci by phenotype after merging
#Merged across all annotations
VF=E4
while read pheno; do
  for dummy in 1; do
    echo ${pheno}
    for CNV in DEL DUP; do
      fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )
#Specific to cancer (i.e. not observed in any germline group)
while read pheno; do
  for CNV in DEL DUP; do
    fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed | \
    cut -f4
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk '{ if ($2=="GERM") print $1 }' ) | sort | uniq | \
fgrep -wvf - \
${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.all_classes.bed | wc -l
#Count of del-only, del+dup, and dup-only
VF=E4
while read pheno; do
  for dummy in 1; do
    echo ${pheno}
    #DEL-only
    fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.all_classes.bed | \
    cut -f4 | fgrep -wvf - \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.all_classes.bed | wc -l
    #DEL+DUP
    fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.all_classes.bed | \
    cut -f4 | fgrep -wf - \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.all_classes.bed | wc -l
    #DUP-only
    fgrep -v "#" ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.all_classes.bed | \
    cut -f4 | fgrep -wvf - \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.all_classes.bed | wc -l
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )







