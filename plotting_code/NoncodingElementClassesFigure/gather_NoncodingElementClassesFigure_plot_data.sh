#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure for noncoding element class-based rCNV burden tests

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure ]; then
  rm -rf ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure
fi
mkdir ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure

#####Get counts of del-only, dup-only, and del+dup sites per 
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
          cut -f1 | fgrep -v CTRL ) > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/significant_loci_per_pheno.del_dup_both.counts.txt

#####Get list of element sizes significant in GERM, NEURO, SOMA, and CNCR
VF=E4
if ! [ -e ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros ]; then
  mkdir ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros
fi
#GERM
while read pheno; do
  for CNV in DEL DUP; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( sed -n '3,24p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
sort | uniq | fgrep -wf - \
${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed | \
awk '{ print $3-$2 }' > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros/GERM.element_sizes.txt
#NEURO
while read pheno; do
  for CNV in DEL DUP; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( sed -n '4,13p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
sort | uniq | fgrep -wf - \
${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed | \
awk '{ print $3-$2 }' > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros/NEURO.element_sizes.txt
#SOMA
while read pheno; do
  for CNV in DEL DUP; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( sed -n '15,24p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
sort | uniq | fgrep -wf - \
${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed | \
awk '{ print $3-$2 }' > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros/SOMA.element_sizes.txt
#CNCR
while read pheno; do
  for CNV in DEL DUP; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( sed -n '25,37p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
sort | uniq | fgrep -wf - \
${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed | \
awk '{ print $3-$2 }' > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_size_distros/CNCR.element_sizes.txt

#####Get jaccard index matrix across organ systems per phenotype
if ! [ -e ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/tissue_jaccards ]; then
  mkdir ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/tissue_jaccards
fi
while read pheno; do
  while read annoA; do
    while read annoB; do
      bedtools jaccard -g /data/talkowski/rlc47/src/GRCh37.genome \
      -a <( cat ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.${annoA}.bed \
                ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.${annoA}.bed | \
            sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | uniq ) \
      -b <( cat ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.${annoB}.bed \
                ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.${annoB}.bed | \
            sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | uniq ) | sed '1d' | awk '{ print $3 }'
    done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
              sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - ) | paste -s
  done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
            sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - ) > \
  ${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/tissue_jaccards/${pheno}.jaccard_matrix.txt
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Get count of significant loci across organ systems per phenotype
VF=E4
while read pheno; do
  while read annoA; do
    cat ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DEL_${VF}.final_merged_loci.${annoA}.bed \
              ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_DUP_${VF}.final_merged_loci.${annoA}.bed | \
    cut -f4 | sort | uniq | wc -l
  done < <( cut -f1 ${WRKDIR}/bin/rCNVmap/misc/OrganGroup_Consolidation_NoncodingAnnotation_Linkers.list | \
    sort | uniq | cat <( echo -e "all_classes\ntissue_agnostic" ) - ) | paste -s 
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL ) > \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/element_count_pheno_by_tissue.matrix.txt

######Get fraction of elements overlapping various classes of noncoding elements
#Create master CNCR and GERM element lists for DEL and DUP
for CNV in DEL DUP; do
  #GERM
  while read pheno; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done < <( sed -n '3,24p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
  sort | uniq | fgrep -wf - ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed > \
  ${TMPDIR}/GERM.${CNV}.all_elements.bed
  #CNCR
  while read pheno; do
    cut -f4 ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done < <( sed -n '25,37p' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 ) | \
  sort | uniq | fgrep -wf - ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed > \
  ${TMPDIR}/CNCR.${CNV}.all_elements.bed
done

#####Create master matrix of all final merged loci vs all noncoding annotation tracks
VF=E4
annoSet=all_classes
#Initialize matrix
cat ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/${annoSet}_haplosuffDELnoncodingDUP_${VF}.signif_loci.merged.filtered.bed > \
${TMPDIR}/master_noncoding_matrix.build.tmp
#Iterate over matrix and annotate if significant in DEL/DUP per phenotype
for CNV in DEL DUP; do
  while read pheno; do
    bedtools intersect -c \
    -a ${TMPDIR}/master_noncoding_matrix.build.tmp \
    -b ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed > \
    ${TMPDIR}/master_noncoding_matrix.build.tmp2
    mv ${TMPDIR}/master_noncoding_matrix.build.tmp2 \
    ${TMPDIR}/master_noncoding_matrix.build.tmp
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )
done
#Iterate over annotations and add number of hits
while read anno; do
  echo ${anno}
  #Get annotation path
  annopath=$( fgrep -w ${anno} ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.list | cut -f2 )
  bedtools intersect -c -a ${TMPDIR}/master_noncoding_matrix.build.tmp \
  -b <( awk -v OFS="\t" '{ if ($1>0 && $3>$2 && $2>0) print $1, $2, $3 }' ${annopath} ) > \
  ${TMPDIR}/master_noncoding_matrix.build.tmp2
  mv ${TMPDIR}/master_noncoding_matrix.build.tmp2 \
  ${TMPDIR}/master_noncoding_matrix.build.tmp
done < <( cat ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list <( echo -e "" ) )
#Add column names
for dummy in 1; do
  echo -e "#chr\tstart\tend\tlocus_ID"
  for CNV in DEL DUP; do
    while read pheno; do
      echo -e "${pheno}_${CNV}_significant"
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )
  done
  cat ${WRKDIR}/bin/rCNVmap/misc/master_noncoding_annotations.alternative_sort.list
done | paste -s > ${WRKDIR}/analysis/perAnno_burden/master_noncoding_loci.annotated_matrix.bed
cat ${TMPDIR}/master_noncoding_matrix.build.tmp >> \
${WRKDIR}/analysis/perAnno_burden/master_noncoding_loci.annotated_matrix.bed
gzip -f ${WRKDIR}/analysis/perAnno_burden/master_noncoding_loci.annotated_matrix.bed
cp ${WRKDIR}/analysis/perAnno_burden/master_noncoding_loci.annotated_matrix.bed.gz \
${WRKDIR}/data/plot_data/NoncodingElementClassesFigure/










