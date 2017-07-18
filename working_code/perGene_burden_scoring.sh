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
if [ -e ${WRKDIR}/data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/perGene_burden
fi
mkdir ${WRKDIR}/data/perGene_burden

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
        while read pheno; do
          cat ${WRKDIR}/analysis/perGene_burden/signif_genes/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.genes.list | \
          sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
                  fgrep -v CTRL | cut -f1 )
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
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) \
          <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) | \
          sed 's/_/\-/g' | sort | uniq > \
          ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
        done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
                  fgrep -v CTRL | cut -f1 )
        #Cases - merged
        # for group in GERM NEURO NDD PSYCH SOMA CNCR; do
        #   fgrep -wvf \
        #   <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) \
        #   <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) | \
        #   sed 's/_/\-/g' | sort | uniq > \
        #   ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
        # done
        # #Controls
        # for group in GERM NEURO NDD PSYCH SOMA CNCR; do
        #   sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_${group}_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list
        # done | sort | uniq | fgrep -wvf - \
        # <( sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.union.genes.list ) > \
        # ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/all_CTRL_groups_${CNV}_${VF}_${context}.geneScore_${sig}_sig.unique.genes.list
      done
    done
  done
done









