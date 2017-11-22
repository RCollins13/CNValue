#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run all rCNV burden scoring per gene

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
# if [ -e ${WRKDIR}/data/perGene_burden ]; then
#   rm -rf ${WRKDIR}/data/perGene_burden
# fi
# mkdir ${WRKDIR}/data/perGene_burden

#####Create subdirectories
while read pheno; do
  if [ -e ${WRKDIR}/data/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perGene_burden/${pheno}
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Submit burden data collection for all phenotypes
#Germline
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      #Exonic
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_exonic \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_exonic.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
      #Wholegene
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_wholegene \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh -W \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_wholegene.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )
#Final analysis: subset of all combinations, but exclude CNVs with >200kb overlap with any large segment
#Prep CNVs:
for pheno in CTRL GERM NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    for VF in E2 E4; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz | head -n1 > \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed
      bedtools coverage \
      -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_E2_all.signif.bed \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all.bed.gz | \
      awk -v OFS="\t" '{ if ($(NF-2)<200000) print $1, $2, $3, $4, $5, $6, $7 }' | \
      sort -Vk1,1 -k2,2n -k3,3n >> \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed
      gzip -f ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed
    done
  done
done
#Launch jobs
for pheno in GERM NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    for VF in E2 E4; do
      #Exonic
      bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_exonic \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_exonic.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
      #Wholegene
      bsub -q normal -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_wholegene \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh -W \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_wholegene.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.all_largeSegmentExcluded.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
    done
  done
done


#####Copy test datasets to plotting data directory
if [ -e ${WRKDIR}/data/plot_data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/plot_data/perGene_burden
fi
mkdir ${WRKDIR}/data/plot_data/perGene_burden
for pheno in GERM CNCR; do
  for CNV in CNV DEL DUP; do
    for VF in E2 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          cp ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          ${WRKDIR}/data/plot_data/perGene_burden/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt
        fi
      done
    done
  done
done

#####Copy all datasets to plotting data directory
if [ -e ${WRKDIR}/data/plot_data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/plot_data/perGene_burden
fi
mkdir ${WRKDIR}/data/plot_data/perGene_burden
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          cp ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          ${WRKDIR}/data/plot_data/perGene_burden/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          cut -f1 | fgrep -v CTRL )

#####Run geneScore model for all combinations
#Create master output directory
if [ -e ${WRKDIR}/analysis/perGene_burden ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden
fi
mkdir ${WRKDIR}/analysis/perGene_burden
#Run model
while read pheno; do
  #Create output directory for results
  if [ -e ${WRKDIR}/analysis/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perGene_burden/${pheno}
  #Get number of subjects in group
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
           ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${context}_geneScoreModel \
          "${WRKDIR}/bin/rCNVmap/bin/run_geneScore_model.R \
          -o ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
          ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          38628 ${nCASE}"
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )
#Final analysis: subset of all combinations
for pheno in GERM NEURO NDD PSYCH SOMA; do
  nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
         ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
  for CNV in DEL DUP; do
    for VF in E2 E4; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${context}_geneScoreModel \
          "${WRKDIR}/bin/rCNVmap/bin/run_geneScore_model.R \
          -o ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt \
          ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt \
          38628 ${nCASE}"
        fi
      done
    done
  done
done


#####Subselect significant genes per combination
#Create master output directory
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes
#Per disease phenotype
while read pheno; do
  #Create output directory for results
  if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          #Nominal
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($43<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_nominally_sig.genes.list
          #FDR
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($44<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_FDR_sig.genes.list
          #Bonferroni
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($45<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_Bonferroni_sig.genes.list
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )
#Final analysis: subset of all combinations
for pheno in GERM NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    for VF in E2 E4; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          #Nominal
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($43<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_nominally_sig.genes.list
          #FDR
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($44<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_FDR_sig.genes.list
          #Bonferroni
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($45<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_Bonferroni_sig.genes.list
        fi
      done
    done
  done
done
#Create output directory for control results
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL
#Genes per control comparison
while read pheno; do  
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for context in exonic wholegene; do
        if [ -e ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_data.txt ]; then
          #Nominal
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($46<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL/CTRL_vs_${pheno}_${CNV}_${VF}_${context}.geneScore_nominally_sig.genes.list
          #FDR
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($47<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL/CTRL_vs_${pheno}_${CNV}_${VF}_${context}.geneScore_FDR_sig.genes.list
          #Bonferroni
          fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt | \
          awk -v FS="\t" '{ if ($48<=0.05) print $1 }' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL/CTRL_vs_${pheno}_${CNV}_${VF}_${context}.geneScore_Bonferroni_sig.genes.list
        fi
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )

#####Get summary table of significant gene counts
VF=E4
CNV=DEL
context=exonic
while read pheno; do
  for dummy in 1; do
    echo ${pheno}
    for sig in nominally FDR Bonferroni; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list | wc -l
    done
    for sig in nominally FDR Bonferroni; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL/CTRL_vs_${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )

#####Get union & unique lists of genes per
#Prep output directory
if [ -e ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/ ]; then
  rm -rf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/
fi
mkdir ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/
#Iterate over all CNV filters
for CNV in CNV DEL DUP; do
  echo ${CNV}
  for VF in E2 E3 E4 N1; do
    echo ${VF}
    for context in exonic wholegene; do
      echo ${context}
      for sig in nominally FDR Bonferroni; do
        echo ${sig}
        #Each original pheno group
        # while read pheno; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list | \
        #   sort | uniq > \
        #   ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
        #           fgrep -v CTRL | cut -f1 )
        # #GERM
        # while read pheno; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done < <( awk '{ if ($2=="GERM") print $1 }' \
        #           ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
        # sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_GERM_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #NEURO
        # for pheno in NEURO NDD ASD DD ID PSYCH SCZ BEHAV SEIZ HYPO; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done | sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_NEURO_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #NDD
        # for pheno in NDD ASD DD ID; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done | sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_NDD_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #PSYCH
        # for pheno in PSYCH SCZ BEHAV; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done | sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_PSYCH_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #SOMA
        # for pheno in SOMA DRU GRO INT SKEL HEAD MSC CARD EE EMI; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done | sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_SOMA_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #CNCR
        # while read pheno; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done < <( awk '{ if ($2=="CNCR") print $1 }' \
        #           ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list ) | \
        # sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CNCR_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # #CTRL
        # while read pheno; do
        #   cat ${WRKDIR}/analysis/perGene_burden/signif_genes/CTRL/CTRL_vs_${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list
        # done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
        #           fgrep -v CTRL | cut -f1 ) | \
        # sort | uniq > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        ##############
        # Unique lists
        ##############
        #Each original pheno group
        while read pheno; do
          fgrep -wvf \
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_*_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) \
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) | \
          sed 's/_/\-/g' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
        done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
                  fgrep -v CTRL | cut -f1 )
        # Cases - merged
        for group in GERM NEURO NDD PSYCH SOMA CNCR; do
          fgrep -wvf \
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_*_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) \
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) | \
          sed 's/_/\-/g' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
        done
        #Controls
        for group in GERM NEURO NDD PSYCH SOMA CNCR; do
          sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_*_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        done | sort | uniq | fgrep -wvf - \
        <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) > \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
      done
    done
  done
done
#Final analysis: subset of all combinations
#Require E4-sig genes also to be nominally sig via Fisher's exact at E2
#Require <50% gene body overlap versus control probe deserts
for pheno in GERM NEURO NDD PSYCH SOMA; do
  for CNV in DEL DUP; do
    for VF in E4; do
      for context in exonic wholegene; do
        #First must pass E2 basic Fisher's exact, then must pass E4 FDR, then filter against control probe deserts
        fgrep -v "#" ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_E2_${context}.geneScore_stats.txt | \
        awk -v FS="\t" '{ if ($29<=0.05) print $1 }' | sort | uniq | sed 's/\-/_/g' | \
        fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_FDR_sig.genes.list ) | \
        fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.all.bed ) | \
        bedtools coverage -b - -a ${WRKDIR}/data/misc/CTRL_probeDeserts.bed | \
        awk -v OFS="\t" '{ if ($NF<0.5) print $4 }' | sort | uniq > \
        ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list
      done
    done
  done
done
#Count number of genes per phenotype
for pheno in GERM NEURO NDD PSYCH SOMA; do
  echo ${pheno}
  for CNV in DEL DUP; do
    for VF in E4; do
      for context in exonic; do
        cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list | wc -l
      done
    done
  done
done | paste - - -

#Sanity check: get fraction of constrained genes per group for final analysis set
VF=E4; context=exonic
for CNV in DEL DUP; do
  for context in exonic wholegene; do
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      for wrapper in 1; do
        echo "${pheno}_${CNV}_${context}"
        cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list | wc -l
        for geneSet in ExAC_constrained ExAC_missense_constrained DDD_2017 DDG2P_AnyConf_Dominant_LOF DDG2P_AnyConf_Dominant_GOF DDG2P_AnyConf_Dominant_Unknown; do
          total=$( cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list | wc -l )
          hit=$( fgrep -wf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list \
                 ${WRKDIR}/data/master_annotations/genelists/${geneSet}.genes.list | wc -l )
          echo "" | awk -v hit=${hit} -v total=${total} '{ print hit/total }'
        done
      done | paste -s
    done
  done 
done
for geneSet in ExAC_constrained ExAC_missense_constrained DDD_2017 DDG2P_AnyConf_Dominant_LOF DDG2P_AnyConf_Dominant_GOF DDG2P_AnyConf_Dominant_Unknown; do
  total=$( cat ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | wc -l )
  hit=$( fgrep -wf ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
         ${WRKDIR}/data/master_annotations/genelists/${geneSet}.genes.list | wc -l )
  echo "" | awk -v hit=${hit} -v total=${total} '{ print hit/total }'
done | paste -s

# #####Get summary table of significant gene counts - genes unique to disease or control
# VF=E4
# context=exonic
# for CNV in CNV DEL DUP; do
#   echo -e "\n\n${CNV}\n----"
#   for dummy in 1; do
#     #Control
#     for dummy in 1; do
#       echo "CTRL"
#       for sig in nominally FDR Bonferroni; do
#         cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
#       done
#     done | paste -s
#     #Affected pheno groups
#     while read pheno; do
#       for dummy in 1; do
#         echo ${pheno}
#         for sig in nominally FDR Bonferroni; do
#           cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
#         done
#       done | paste -s
#     done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
#               fgrep -v CTRL | cut -f1 )
#     #Affected - unions
#     for group in GERM NEURO NDD PSYCH SOMA CNCR; do
#       for dummy in 1; do
#         echo "ALL_${group}_UNION"
#         for sig in nominally FDR Bonferroni; do
#           cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list | wc -l
#         done
#       done | paste -s
#     done
#   done
# done

# #####Get count of germline pheno groups significant per gene (>0)
# VF=E4
# context=exonic
# sig=Bonferroni
# while read pheno; do
#   for CNV in DEL DUP; do
#     cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
#   done | sort | uniq
# done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
#           awk '{ if ($2=="GERM") print $1 }' ) | \
# sort | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | sort -nrk2,2


#####Final analysis: get list of all del-significant, dup-significant, and any-significant E4 exonic genes in top 5 phenos
#MANUALLY REMOVE GENES WHERE REGULATORY EFFECT IS MORE LIKELY
#SEE perAnno_burden_scoring.sh FOR DETAILS
#GENES TO EXCLUDE: CDH19, ERICH1, GNAT3, GPC5, GTDC1, LRRFIP2, OVGP1, SEL1L2, ZNF257, ZNF676
VF=E4
for context in exonic; do
  for CNV in DEL DUP; do
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      echo -e "CDH19\nERICH1\nGNAT3\nGPC5\nGTDC1\nLRRFIP2\nOVGP1\nSEL1L2\nZNF257\nZNF676" | \
      fgrep -wvf - ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list > \
      ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list2
      mv ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list2 \
      ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list
    done
  done
done

#DEL OR DUP SEPARATELY
VF=E4
for context in exonic; do
  for CNV in DEL DUP; do
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list
    done | sort | uniq > ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list
  done
done
#DEL AND DUP UNION
for context in exonic; do
  for CNV in DEL DUP; do
    cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list
  done | sort | uniq > ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_DELDUPUnion_${VF}_${context}.geneScore_FINAL_sig.genes.list
done

#####Get genes with no prior germline disease association
while read gene; do
  hits=$( fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/* | \
          fgrep -v GTEx | fgrep -v MASTER | fgrep -v GO_ | fgrep -v quintile | \
          fgrep -v decile | fgrep -v targets | fgrep -v cancer | fgrep -v Cancer | \
          fgrep -v Gencode | fgrep -v ExAC | fgrep -v Constrained | fgrep -v RVIS | \
          fgrep -v essential | fgrep -v DNA_repair | fgrep -v GWAS_nearest | \
          fgrep -v Olfactory | fgrep -v COSMIC | fgrep -v Kinases | fgrep -v GPCRs | sed '/^$/d' | wc -l )
  echo -e "${gene}\t${hits}"
done < ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_DELDUPUnion_${VF}_${context}.geneScore_FINAL_sig.genes.list
#Read out other gene lists for constrained genes with no disease associations
while read gene; do
  echo -e "\n\n${gene}"
  fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/* | cut -f1 -d\: | \
  awk -v FS="/" '{ print $NF }' | sed 's/\.genes\.list//g' | fgrep -v Gencode | \
  fgrep -v GTEx
done < <( echo -e "ACTR3B\nADCY9\nANKFY1\nASB7\nBPTF\nC16orf72\nCLUH\nCOLEC12\nFAM196A\nFCGR1B\nFNTA\nGMDS\nGPR89A\nLRRK1\nMPRIP\nNCS1\nOR4F17\nOR4F4\nOTUD7A\nRAP1GAP2\nRAP1GDS1\nRMND5A\nRPSAP58\nTAOK3\nUBAP1\nUSP34\nZNF195" )




