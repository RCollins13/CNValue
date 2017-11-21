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
#Association data
VF=E2
filt=all
for pheno in GERM NEURO NDD PSYCH SOMA; do
  paste \
  <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_DEL_${VF}_${filt}.TBRden_results.bed.gz | \
     sed '1d' | awk -v OFS="\t" '{ print $1, $2, $3, $(NF-1) }' ) \
  <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_DUP_${VF}_${filt}.TBRden_results.bed.gz | \
     sed '1d' | awk -v OFS="\t" '{ print $(NF-1) }' ) > \
  ${WRKDIR}/data/plot_data/figure2/${pheno}.${VF}.${filt}.manhattan_pvals.bed
done
#Novel hits
bedtools coverage \
-a ${WRKDIR}/data/misc/known_CNV_segments/all.merged.bed \
-b ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed | \
awk '{ if ($NF==0) print $0 }' | cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/plot_data/figure2/novel_loci.bed


#####Copy positive control dataset overlaps
cp ${WRKDIR}/analysis/large_CNV_segments/known_locus_overlap.txt \
${WRKDIR}/data/plot_data/figure2/

# #####Gather data for constrained gene enrichment analysis
# cp -r ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments \
# ${WRKDIR}/data/plot_data/figure2/

# #####Gather germline prevalance data
# cp -r ${WRKDIR}/analysis/large_CNV_segments/prevalance \
# ${WRKDIR}/data/plot_data/figure2/