#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for all plots used in figure 2

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/figure2 ]; then
  rm -rf ${WRKDIR}/data/plot_data/figure2
fi
mkdir ${WRKDIR}/data/plot_data/figure2

#####Copy association statistics for CNV loci for E2 all
VF=E2
filt=all
for CNV in DEL DUP; do
  for stat in pVal OR; do
    cp ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_${stat}.bed \
    ${WRKDIR}/data/plot_data/figure2/
  done
done

#####Gather data for manhattan plots
VF=E2
filt=all
for pheno in GERM NEURO NDD PSYCH SOMA; do
  paste \
  <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_DEL_${VF}_${filt}.TBRden_results.bed.gz | \
     sed '1d' | awk -v OFS="\t" '{ print $1, $2, $3, $NF }' ) \
  <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_DUP_${VF}_${filt}.TBRden_results.bed.gz | \
     sed '1d' | awk -v OFS="\t" '{ print $NF }' ) > \
  ${WRKDIR}/data/plot_data/figure2/${pheno}.${VF}.${filt}.manhattan_pvals.bed
done

# #####Gather data for constrained gene enrichment analysis
# cp -r ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments \
# ${WRKDIR}/data/plot_data/figure2/

# #####Gather germline prevalance data
# cp -r ${WRKDIR}/analysis/large_CNV_segments/prevalance \
# ${WRKDIR}/data/plot_data/figure2/