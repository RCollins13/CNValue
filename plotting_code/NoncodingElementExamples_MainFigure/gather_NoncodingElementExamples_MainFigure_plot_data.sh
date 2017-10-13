#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure for specific/individual noncoding elements

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if ! [ -e ${WRKDIR}/data/plot_data/IndividualNoncodingElementsFigure ]; then
  mkdir ${WRKDIR}/data/plot_data/IndividualNoncodingElementsFigure
fi

####Get list of elements ±500kb or E-P connected to (but not overlapping) high-conf disease genes
####NDD
#Get list of genes
cat ${WRKDIR}/data/master_annotations/genelists/extTADA*.genes.list \
${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list | sort | uniq | \
sed -e 's/\-/_/g' -e 's/\./_/g' > ${TMPDIR}/NDD_union.genes.list
#Get 750 flanks and any relevant E-P connections
flank=750000
cat \
<( fgrep -wf ${TMPDIR}/NDD_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) | \
   awk -v OFS="\t" -v flank=${flank} '{ print $1, $2-flank, $2, $4"\n"$1, $3, $3+flank, $4 }' | \
   awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3, $4, "flank" }' ) \
<( fgrep -wf ${TMPDIR}/NDD_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_4.unique_EP_pairs.bed ) | \
   awk -v OFS="\t" '{ print $1, $2, $3, $4, "EP" }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 > \
${TMPDIR}/NDD_union.loci.bed
#Intersect versus significant NDD elements
VF=E4
while read pheno; do
  for CNV in DEL DUP; do
    awk -v OFS="\t" -v pheno=${pheno} -v CNV=${CNV} '{ print $1, $2, $3, $4, pheno"_"CNV }' \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk '{ if ($NF>=10000 && $2=="GERM") print $1 }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | \
bedtools merge -c 4 -o distinct -i - | bedtools intersect -wa -wb -a - \
-b ${TMPDIR}/NDD_union.loci.bed | bedtools intersect -v -wa -a - \
-b <( fgrep -wf ${TMPDIR}/NDD_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) ) > \
${TMPDIR}/NDD_loci_signif.bed
#####CNCR
#Get list of genes
cat ${WRKDIR}/data/master_annotations/genelists/COSMIC*.genes.list | sort | uniq | \
sed -e 's/\-/_/g' -e 's/\./_/g' > ${TMPDIR}/CNCR_union.genes.list
#Get 750kb flanks and any relevant E-P connections
flank=750000
cat \
<( fgrep -wf ${TMPDIR}/CNCR_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) | \
   awk -v OFS="\t" -v flank=${flank} '{ print $1, $2-flank, $2, $4"\n"$1, $3, $3+flank, $4 }' | \
   awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3, $4, "flank" }' ) \
<( fgrep -wf ${TMPDIR}/CNCR_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_4.unique_EP_pairs.bed ) | \
   awk -v OFS="\t" '{ print $1, $2, $3, $4, "EP" }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 > \
${TMPDIR}/CNCR_union.loci.bed
#Intersect versus significant CNCR elements
VF=E4
while read pheno; do
  for CNV in DEL DUP; do
    awk -v OFS="\t" -v pheno=${pheno} -v CNV=${CNV} '{ print $1, $2, $3, $4, pheno"_"CNV }' \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk '{ if ($NF>=10000 && $2=="CNCR") print $1 }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | \
bedtools merge -c 4,5 -o distinct -i - | bedtools intersect -wa -wb -a - \
-b ${TMPDIR}/CNCR_union.loci.bed | bedtools intersect -v -wa -a - \
-b <( fgrep -wf ${TMPDIR}/CNCR_union.genes.list \
   <( sed -e 's/\-/_/g' -e 's/\./_/g' \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) ) > \
${TMPDIR}/CNCR_loci_signif.bed

#####Get number of genes shared between GERM & CNCR from proximity analysis
fgrep -wf <( awk '{ print $(NF-1) }' ${TMPDIR}/CNCR_loci_signif.bed ) \
<( awk '{ print $(NF-1) }' ${TMPDIR}/NDD_loci_signif.bed )

#####Get significant elements with ≥3 gene contacts
####NDD
VF=E4
while read pheno; do
  for CNV in DEL DUP; do
    awk -v OFS="\t" -v pheno=${pheno} -v CNV=${CNV} '{ print $1, $2, $3, $4, pheno"_"CNV }' \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk '{ if ($NF>=10000 && $2=="GERM") print $1 }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | \
bedtools merge -c 4 -o distinct -i - | bedtools intersect -c -a - \
-b ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_1.unique_EP_pairs.bed | \
awk '{ if ($NF>2) print $0 }' > ${TMPDIR}/NDD_multi_signif.loci.bed
####CNCR
VF=E4
while read pheno; do
  for CNV in DEL DUP; do
    awk -v OFS="\t" -v pheno=${pheno} -v CNV=${CNV} '{ print $1, $2, $3, $4, pheno"_"CNV }' \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.all_classes.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk '{ if ($NF>=10000 && $2=="CNCR") print $1 }' ) | \
sort -Vk1,1 -k2,2n -k3,3n -Vk4,4 | \
bedtools merge -c 4 -o distinct -i - | bedtools intersect -c -a - \
-b ${WRKDIR}/data/misc/pcHiC_contacts/pcHiC_contacts.formatted_wGenes.min_1.unique_EP_pairs.bed | \
awk '{ if ($NF>2) print $0 }' > ${TMPDIR}/CNCR_multi_signif.loci.bed

#####Get number of significant elements in both GERM & CNCR
cut -f4 ${TMPDIR}/NDD_multi_signif.loci.bed | fgrep -wf - ${TMPDIR}/CNCR_multi_signif.loci.bed







