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

















