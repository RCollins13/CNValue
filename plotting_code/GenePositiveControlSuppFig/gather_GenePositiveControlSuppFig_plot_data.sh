#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for plots of positive control gene sets for per-gene burden tests

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/GenePositiveControlSuppFig ]; then
  rm -rf ${WRKDIR}/data/plot_data/GenePositiveControlSuppFig
fi
mkdir ${WRKDIR}/data/plot_data/GenePositiveControlSuppFig

#####Step 1: define universe of genes eligible for testing
#Filters: protein-coding, autosomal, <=50% coverage by probe desert
bedtools coverage -a ${WRKDIR}/data/misc/CTRL_probeDeserts.bed \
-b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
awk -v OFS="\t" -v OFS="\t" '{ if ($NF<=0.5) print $1, $4 }' | grep -e '^[0-9]' | \
cut -f2 | sed 's/\-/_/g' | sort | uniq > ${TMPDIR}/gene_universe.txt

#####Step 2: iterate over DEL/DUP, count observed overlap and expected overlap based on fraction of genes in list
VF=E4; context=exonic
for CNV in DEL DUP; do
  for list in ExAC_constrained ExAC_extremely_constrained ExAC_missense_constrained ExAC_haplosufficient \
              ClinGen_haploinsufficient_low_confidence ClinVar_disease_associated \
              DDD_2017 extTADA_DD ASD_TADA_q0.3 extTADA_ASD extTADA_ID \
              Coe2017_ASD_ID_denovolyzeR_LGD Coe2017_ASD_ID_denovolyzeR_MIS \
              Autosomal_dominant_disease Autosomal_recessive_disease \
              DDG2P_HighConf_Dominant_LOF DDG2P_AnyConf_Dominant_LOF DDG2P_HighConf_Recessive_LOF DDG2P_AnyConf_Recessive_LOF \
              DDG2P_HighConf_Dominant_GOF DDG2P_AnyConf_Dominant_GOF \
              DDG2P_HighConf_Dominant_Unknown DDG2P_AnyConf_Dominant_Unknown DDG2P_HighConf_Recessive_Unknown DDG2P_AnyConf_Recessive_Unknown; do
    for wrapper in 1; do
      echo -e "${CNV}\t${list}"
      #Sig genes (all)
      cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list | wc -l
      #Sig genes (in set)
      if [ -e ${WRKDIR}/data/master_annotations/genelists/${list}.genes.list ]; then
        sed 's/\-/_/g' ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list | \
        fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/${list}.genes.list | \
                       fgrep -wf ${TMPDIR}/gene_universe.txt ) | wc -l
      else
        echo NA
      fi
      #Universe genes (all)
      cat ${TMPDIR}/gene_universe.txt | wc -l
      #Universe genes (in set)
      if [ -e ${WRKDIR}/data/master_annotations/genelists/${list}.genes.list ]; then
        sed 's/\-/_/g' ${TMPDIR}/gene_universe.txt | \
        fgrep -wf - <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/genelists/${list}.genes.list | \
                       fgrep -wf ${TMPDIR}/gene_universe.txt ) | wc -l
      else
        echo NA
      fi
    done | paste -s
  done 
done 