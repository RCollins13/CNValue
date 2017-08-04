#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure 4

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/figure4 ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure4
fi
mkdir ${WRKDIR}/data/plot_data/figure4

#####Copy geneScore results for GERM/NEURO/NDD/PSYCH/SOMA/CNCR, exonic/wholegene, E3/E4/N1, CNV/DEL/DUP
#Initialize directory
if [ -e ${WRKDIR}/data/plot_data/figure4/geneScore_results ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure4/geneScore_results
fi
mkdir ${WRKDIR}/data/plot_data/figure4/geneScore_results
#Copy all files
for pheno in GERM NEURO NDD ASD PSYCH SOMA HEAD CNCR; do
  for CNV in CNV DEL DUP; do
    for VF in E3 E4 N1; do
      for context in exonic wholegene; do
        cp ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
        ${WRKDIR}/data/plot_data/figure4/geneScore_results/
      done
    done
  done
done

#####Cut ExAC lof obs/exp, Z-scores, and pLI for correlation vs CNV Z-score
#Constraint
sed '1d' ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | \
awk -v OFS="\t" '{ if ($16>0) print $2, $13/$16, $19, $20 }' > \
${WRKDIR}/data/plot_data/figure4/ExAC_LoF_constraint.txt
#RVIS (0.01%)
sed '1d' ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | \
awk -v OFS="\t" '{ print $1, $(NF-1), $NF }' > \
${WRKDIR}/data/plot_data/figure4/ExAC_RVIS.txt

#####Gather data per phenotype for gene set enrichment plots
#Make universal list of genes tested
fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/GERM/GERM_CNV_E2_exonic.geneScore_stats.txt | \
cut -f1 > ${TMPDIR}/all_tested_genes.list
#Create directories
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/
#Submit data collection - merged master pheno groups
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/geneSet_enrichment_comparison_sets.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/${CNV}_${VF}_${context}_${sig}.comparisons.txt"
      done
    done
  done
done
#All original pheno groups, all gene sets
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons_allPhenos_allGenes -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.all_phenos.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/${CNV}_${VF}_${context}_${sig}.comparisons.txt"
      done
    done
  done
done
#All original pheno groups, all gene sets, exclude constrained genes
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for context in exonic wholegene; do
      for sig in nominally FDR Bonferroni; do
        bsub -q short -sla miket_sc -J ${CNV}_${VF}_${context}_${sig}_collectComparisons_allPhenos_allGenes_noConstrained -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_geneScore_signif_overlaps_vs_other_gene_sets.all_phenos.wExclusion.sh \
        ${CNV} ${VF} ${context} ${sig} \
        ${WRKDIR}/bin/rCNVmap/misc/master_gene_sets.sorted.list \
        ${TMPDIR}/all_tested_genes.list \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/${CNV}_${VF}_${context}_${sig}.comparisons.txt \
        ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list"
      done
    done
  done
done
#Copy to plot data directory
if [ -e ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons
fi
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/subset/
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
mkdir ${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets_noConstrained/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/subset/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/
cp ${WRKDIR}/analysis/perGene_burden/signif_genes/geneset_comparisons_allPhenos_noConstrained/* \
${WRKDIR}/data/plot_data/signif_genes_geneset_comparisons/allPhenos_allSets_noConstrained/

#####Copy significant genes to plotting data directory
#Make directory
if [ -e ${WRKDIR}/data/plot_data/signif_genes_unique ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_genes_unique
fi
mkdir ${WRKDIR}/data/plot_data/signif_genes_unique
cp /data/talkowski/Samples/rCNVmap/analysis/perGene_burden/signif_genes/merged/*unique.genes.list \
${WRKDIR}/data/plot_data/signif_genes_unique/

#####Get union set of NDD and CNCR E4 exonic Bonf sig genes
VF="E4"
for CNV in DEL DUP; do
  for pheno in CNCR NDD; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list
  done
done | sort | uniq > \
${TMPDIR}/NDD_CNCR_DEL_DUP_Bonferroni.genes.list
while read gene; do
  for dummy in 1; do
    echo ${gene}
    for pheno in NDD CNCR; do
      for CNV in DEL DUP; do
        fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list | wc -l
      done
    done
  done | paste -s
done < ${TMPDIR}/NDD_CNCR_DEL_DUP_Bonferroni.genes.list > \
${WRKDIR}/data/plot_data/figure4/NDD_CNCR_PPI_CNVmembership.txt
for CNV in CNV DEL DUP; do
  for pheno in CNCR NDD; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list
  done
done | sort | uniq > \
${TMPDIR}/NDD_CNCR_CNV_DEL_DUP_Bonferroni.genes.list
while read gene; do
  for dummy in 1; do
    echo ${gene}
    for pheno in NDD CNCR; do
      for CNV in CNV DEL DUP; do
        fgrep -w ${gene} ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_exonic.geneScore_Bonferroni_sig.unique.genes.list | wc -l
      done
    done
    fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list | wc -l
  done | paste -s
done < ${TMPDIR}/NDD_CNCR_CNV_DEL_DUP_Bonferroni.genes.list > \
${WRKDIR}/data/plot_data/figure4/NDD_CNCR_PPI_CNVmembership_CNV_DEL_DUP.txt

#####Get counts of significant genes (DEL/DUP) per phenotype
VF=E4
context=exonic
sig=Bonferroni
for dummy in 1; do
  #Control
  echo "CTRL"
  #DEL only
  fgrep -wvf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
  ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
  #Both DEL and DUP
  fgrep -wf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
  ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
  #DUP only
  fgrep -wvf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
  ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
  #Affected pheno groups
  while read pheno; do
    echo "${pheno}"
    #DEL only
    fgrep -wvf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
    ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
    #Both DEL and DUP
    fgrep -wf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
    ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
    #DUP only
    fgrep -wvf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list \
    ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
            fgrep -v CTRL | cut -f1 )
done | paste - - - - > \
${WRKDIR}/data/plot_data/figure4/signif_genes_by_pheno.count.txt

#####Get union of all genes signif in any disease group and not controls
VF=E4
context=exonic
sig=Bonferroni
while read pheno; do
  for CNV in DEL DUP; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 ) | \
sort | uniq

#####Get union of all genes signif in any control comparison and not any disease groups
VF=E4
context=exonic
sig=Bonferroni
for CNV in DEL DUP; do
  cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
done | sort | uniq

#####Print list of genes significant in any cancer or any germline group (DEL/DUP exonic E4 only)
VF=E4
context=exonic
sig=Bonferroni
for group in GERM CNCR; do
  while read pheno; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done | sort | uniq
  done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort | uniq -c | awk '{ if ($1>=3) print $2 }' > ${TMPDIR}/${group}_union_genes.list
done
#GERM union vs CNCR union - asterisk if constrained
while read gene; do
  const=$( fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | wc -l )
  if [ ${const} -gt 0 ]; then
    echo "${gene}*"
  else
    echo ${gene}
  fi
done < <( fgrep -wf <( sed 's/\-/_/g' ${TMPDIR}/GERM_union_genes.list ) \
          <( sed 's/\-/_/g' ${TMPDIR}/CNCR_union_genes.list ) | \
          sed 's/\_/\-/g' | sort | uniq )

#####Print list of pleiotropic genes significant in at least N phenotype groups (1 cancer/1 germ min)
VF=E4
context=exonic
sig=Bonferroni
N=10
for group in GERM CNCR; do
  while read pheno; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done | sort | uniq
  done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort | uniq -c | awk -v OFS="\t" -v N=${N} '{ if ($1>=N) print $2, $1 }'
done | fgrep -wf ${TMPDIR}/CNCR_union_genes.list | \
fgrep -wf ${TMPDIR}/GERM_union_genes.list | sort -nrk2,2 | cut -f1 > \
${TMPDIR}/pleiotropicGenes_CNCR_GERM.genes.list
#Write header
for dummy in 1; do
  echo "gene"
  sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v "CTRL" | cut -f1
done | paste -s > \
${WRKDIR}/data/plot_data/figure4/pleiotropicGenes_CNCR_GERM.genes.list
#Write contents
while read gene; do
  for dummy in 1; do
    #asterisk if constrained
    const=$( fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
    <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | wc -l )
    if [ ${const} -gt 0 ]; then
      echo -e "${gene}*"
    else
      echo -e "${gene}"
    fi
    #Print support from each disease
    while read pheno; do
      DEL=$( echo "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
             <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list ) | \
             wc -l )
      DUP=$( echo "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
             <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list ) | \
             wc -l )
      if [ ${DEL} -gt 0 ] && [ ${DUP} -gt 0 ]; then
        echo BOTH
      elif [ ${DEL} -gt 0 ] && [ ${DUP} -eq 0 ]; then
        echo DEL
      elif [ ${DEL} -eq 0 ] && [ ${DUP} -gt 0 ]; then
        echo DUP
      else
        echo NOT
      fi
    done < <( sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v "CTRL" | cut -f1 )
  done | paste -s
done < ${TMPDIR}/pleiotropicGenes_CNCR_GERM.genes.list >> \
${WRKDIR}/data/plot_data/figure4/pleiotropicGenes_CNCR_GERM.genes.list

#####Get number of significant groups per gene for DEL/DUP
VF=E4
context=exonic
sig=Bonferroni
N=6
#Get list of all significant genes in >N pheno groups for DEL/DUP
for CNV in DEL DUP; do
  while read pheno; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
  done < <( sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
    fgrep -v CTRL | cut -f1 ) | sort | uniq -c | awk -v N=${N} '{ if ($1>N) print $2 }' > \
  ${TMPDIR}/N${N}_Sig_${CNV}.genes.list
done
while read gene; do
  #asterisk if constrained
  const=$( fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | wc -l )
  if [ ${const} -gt 0 ]; then
    echo -e "${gene}*"
  else
    echo -e "${gene}"
  fi
  for group in GERM CNCR; do
    for CNV in DEL DUP; do
      while read pheno; do
        echo ${gene} | sed 's/\-/_/g' | fgrep -wf - \
        <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list ) | wc -l
      done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
      awk '{ sum+=$1 }END{ print sum }'
    done
  done
done < <( cat ${TMPDIR}/N${N}_Sig_*.genes.list | sort | uniq ) | paste - - - - - > \
  ${WRKDIR}/data/plot_data/figure4/sigGenes_min${N}_countsPerGene_byGroup.txt  

#####Print lists of genes significant in >0 CNCR or >0 GERM
VF=E4
context=exonic
sig=Bonferroni
for group in GERM CNCR; do
  while read pheno; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done | sort | uniq
  done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort -Vk1,1 | uniq > ${TMPDIR}/${group}_union_genes.list
done

#####Print list of pleiotropic genes significant in at least 17/22 (>75%) germ phenos and NO cancers
VF=E4
context=exonic
sig=Bonferroni
N=6
for group in GERM; do
  while read pheno; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done | sort | uniq
  done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort | uniq -c 
done | fgrep -wvf ${TMPDIR}/CNCR_union_genes.list | \
awk -v OFS="\t" -v N=${N} '{ if ($1>=N) print $2, $1 }' | sort -nrk2,2 > \
${TMPDIR}/${N}min_germ_only_genes.list
#asterisk if constrained
while read gene count; do
  const=$( fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | wc -l )
  if [ ${const} -gt 0 ]; then
    echo -e "${gene}*\t${count}"
  else
    echo -e "${gene}\t${count}"
  fi
done < ${TMPDIR}/${N}min_germ_only_genes.list

#####Print list of genes significant in at least 11/13 (>90%) CNCR phenos and NO germline
VF=E4
context=exonic
sig=Bonferroni
N=11
for group in CNCR; do
  while read pheno; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
    done | sort | uniq
  done < <( awk -v group=${group} '{ if ($2==group) print $1 }' \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
  sort | uniq -c 
done | fgrep -wvf ${TMPDIR}/GERM_union_genes.list | \
awk -v OFS="\t" -v N=${N} '{ if ($1>=N) print $2, $1 }' | sort -nrk2,2 > \
${TMPDIR}/${N}min_cncr_only_genes.list
#asterisk if constrained
while read gene count; do
  const=$( fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
  <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list ) | wc -l )
  if [ ${const} -gt 0 ]; then
    echo -e "${gene}*\t${count}"
  else
    echo -e "${gene}\t${count}"
  fi
done < ${TMPDIR}/${N}min_cncr_only_genes.list

#####Get matrix of associations per top genes
#Write header
for dummy in 1; do
  echo "gene"
  sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v "CTRL" | cut -f1
done | paste -s > ${WRKDIR}/data/plot_data/figure4/topGenes_assocByPheno.txt
#Write contents
while read gene; do
  for dummy in 1; do
    echo -e "${gene}"
    while read pheno; do
      DEL=$( fgrep -w ${gene} \
             ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DEL_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      DUP=$( fgrep -w ${gene} \
             ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l )
      if [ ${DEL} -gt 0 ] && [ ${DUP} -gt 0 ]; then
        echo BOTH
      elif [ ${DEL} -gt 0 ] && [ ${DUP} -eq 0 ]; then
        echo DEL
      elif [ ${DEL} -eq 0 ] && [ ${DUP} -gt 0 ]; then
        echo DUP
      else
        echo NOT
      fi
    done < <( sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v "CTRL" | cut -f1 )
  done | paste -s
done < <( cat ${TMPDIR}/10min_union_genes.list \
              <( echo "SKIP" ) \
              ${TMPDIR}/17min_germ_only_genes.list \
              <( echo "SKIP" ) \
              ${TMPDIR}/11min_cncr_only_genes.list | cut -f1 ) >> \
${WRKDIR}/data/plot_data/figure4/topGenes_assocByPheno.txt

######Process Tanya's GO term enrichment analyses
#Get list of all terms with Bonf-sig enrichment for DEL/DUP
#N=4,807 total terms tested (correct at 0.05/(35*4807) )
for CNV in DEL DUP; do
  while read pheno; do
    fgrep "\"" \
    ${WRKDIR}/analysis/GO_enrichments/perGene_burden/merged/${pheno}_${CNV}_E4_exonic.geneScore_Bonferroni_sig.unique-classic-BP.txt | \
    sed '1d' | sed 's/\"//g' | \
    awk -v pheno=${pheno} -v FS="\t" -v OFS="\t" '{ if ($6<=0.001) print pheno, $0 }'
  done < <( fgrep -v "#" \
            ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
            fgrep -v CTRL | cut -f1 ) | cut -f2 | sort | uniq
done
"GO:0007399"

