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
for CNV in DEL DUP; do
  for list in ExAC_constrained ExAC_extremely_constrained ExAC_missense_constrained ExAC_haplosufficient \
              ClinGen_haploinsufficient_low_confidence ClinVar_disease_associated \
              DDD_2017 extTADA_DD ASD_TADA_q0.3 extTADA_ASD extTADA_ID \
              Coe2017_ASD_ID_denovolyzeR_LGD Coe2017_ASD_ID_denovolyzeR_MIS \
              extTADA_SEIZ extTADA_SCZ \
              Autosomal_dominant_disease Autosomal_recessive_disease \
              DDG2P_HighConf_Dominant_LOF DDG2P_AnyConf_Dominant_LOF DDG2P_HighConf_Recessive_LOF DDG2P_AnyConf_Recessive_LOF \
              DDG2P_HighConf_Dominant_GOF DDG2P_AnyConf_Dominant_GOF DDG2P_HighConf_Recessive_GOF DDG2P_AnyConf_Recessive_GOF; do
  done
done 