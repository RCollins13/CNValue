#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to compile all supplementary tables for rCNV paper

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Create directory if necessary
if ! [ -e ${WRKDIR}/data/plot_data/suppTables ]; then
  mkdir ${WRKDIR}/data/plot_data/suppTables
fi

#####Supp tables 1-2: significant large segments
VF=E2; filt=all
for CNV in DEL DUP; do
  #Print header
  for wrapper in 1; do
    echo -e "chr\tstart\tend\tCTRL_${CNV}_count"
    for pheno in GERM NEURO NDD PSYCH NONN; do
      echo -e "${pheno}_${CNV}_count\t${pheno}_${CNV}_OR\t${pheno}_${CNV}_p"
    done
  done | paste -s > ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt
  #Get ORs and p-values
  bedtools intersect -wa -wb \
  -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR.bed ) | \
  bedtools intersect -wa -wb -a - \
  -b <( sed 's/^chr//g' ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal.bed ) | \
  cut --complement -f4-6,12-14 > \
  ${TMPDIR}/${CNV}_segments_stats.txt
  #Get count of CNVs per pheno
  while read chr start end; do
    for pheno in CTRL GERM NEURO NDD PSYCH SOMA; do
      bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
      cut -f4- | bedtools coverage -b - -a <( echo -e "${chr}\t${start}\t${end}" ) | \
      awk -v overlap=${overlap} '{ if ($(NF-2)>=overlap) print $4 }' | sort | uniq | wc -l
    done | paste -s
  done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed > \
  ${TMPDIR}/${CNV}_segments_CNVcounts.txt
  #Paste the two tables together
  paste ${TMPDIR}/${CNV}_segments_stats.txt ${TMPDIR}/${CNV}_segments_CNVcounts.txt | \
  awk -v OFS="\t" '{ print $1, $2, $3, $14, $15, $4, $9, $16, $5, $10, $17, $6, $11, $18, $7, $12, $19, $8, $13 }' >> \
  ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt
done

#####Supp tables 3-4: significant genes
VF=E4; context=exonic
for CNV in DEL DUP; do
  #Print header
  for wrapper in 1; do
    echo -e "gene\tCTRL_${CNV}_count"
    for pheno in GERM NEURO NDD PSYCH NONN; do
      echo -e "${pheno}_${CNV}_count\t${pheno}_${CNV}_OR\t${pheno}_Zscore\t${pheno}_${CNV}_q"
    done
  done | paste -s > ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt
  #Iterate over genes & get association statistics per phenotype
  while read gene; do
    for wrapper in 1; do
      echo "${gene}"
      for pheno in GERM NEURO NDD PSYCH SOMA; do
        awk -v OFS="\t" -v gene=${gene} '{ if ($1==gene) print $8, $14, $26, $41, $44 }' \
        ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt
      done
    done | paste -s | cut --complement -f7,12,17,22
  done < ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list >> \
  ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt
done

#####Supp tables 5-6: significant regulatory blocks
VF=E4; annoSet=all_classes
for CNV in DEL DUP; do
  #Print header
  for wrapper in 1; do
    echo -e "chr\tstart\tend\tblockID\tsignif_elements\tCTRL_${CNV}_count"
    for pheno in GERM NEURO NDD PSYCH NONN; do
      echo -e "${pheno}_${CNV}_count\t${pheno}_${CNV}_combined_OR\t${pheno}_${CNV}_minimum_p"
    done
  done | paste -s > ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt
  #Get list of all DEL or DUP significant elements
  for pheno in GERM NEURO NDD PSYCH SOMA; do
    cat ${WRKDIR}/analysis/perAnno_burden/signif_elements/all_merged/final_loci/${pheno}/${pheno}_${CNV}_${VF}.final_merged_loci.${annoSet}.bed
  done | sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 | bedtools merge -i - > \
  ${TMPDIR}/${CNV}_sig_elements.bed
  #Get CNV counts, ORs, and p-values
  bedtools intersect -wa -u -b ${TMPDIR}/${CNV}_sig_elements.bed \
  -a <( sed 's/chr//g' ${WRKDIR}/analysis/perAnno_burden/signifRegulatoryBlocks.${CNV}_ORs.bed ) > \
  ${TMPDIR}/${CNV}_sig_elements.stats.bed 
  #Get individual elements to append to tracks
  while read chr start end blockID everythingElse; do
    for wrapper in 1; do
      echo -e "${chr}\t${start}\t${end}\t${blockID}"
      fgrep -w ${blockID} ${WRKDIR}/analysis/perAnno_burden/signifRegulatoryBlocks.final.bed | cut -f5
      echo "${everythingElse}"
    done | paste -s
  done < ${TMPDIR}/${CNV}_sig_elements.stats.bed >> \
  ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt
done







