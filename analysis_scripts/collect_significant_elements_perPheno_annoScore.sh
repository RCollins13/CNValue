#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects significant elements per annotation class

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reads arguments
pheno=$1
CNV=$2
VF=$3
filt=$4
list=$5

#####Iterates over elements list and runs model
while read class path; do
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}
  fi
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}
  fi
  STATS=${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScore_stats.bed.gz
  if [ -e ${STATS} ]; then
    #Nominally significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($44<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.nom_sig_elements.bed
    #BH-corrected significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($45<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.BH_sig_elements.bed
    #Bonferroni-corrected significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($46<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.bonf_sig_elements.bed
    #Gzip all
    gzip -f ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.*_sig_elements.bed*

    #Nominally significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($47<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.nom_sig_elements.bed
    #BH-corrected significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($48<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.BH_sig_elements.bed
    #Bonferroni-corrected significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($49<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.bonf_sig_elements.bed
    #Gzip all
    gzip -f ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.*_sig_elements.bed*
  else
    echo "OUTPUT FOR ${pheno} ${CNV} ${VF} ${filt} vs ${path} IS MISSING. SKIPPING..."
  fi
done < ${list}