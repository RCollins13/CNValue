#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Gets lowest p-values for a list of loci from the CNV hotspot analysis

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
VF=$1
filt=$2
CNV=$3

#####Run
while read chr start end; do
  for dummy in 1; do
    echo -e "chr${chr}\n${start}\n${end}"
    #iterate over phenos
    while read pheno; do
      bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      -b <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_${CNV}_${VF}_${filt}.TBRden_results.bed.gz | \
      sed '1d' ) | awk '{ print $NF }' > ${TMPDIR}/${pheno}_${CNV}_${VF}_${filt}_chr${chr}_${start}_${end}.pvals.tmp
      unset R_HOME
      Rscript -e "cat(paste(format(min(read.table(\"${TMPDIR}/${pheno}_${CNV}_${VF}_${filt}_chr${chr}_${start}_${end}.pvals.tmp\",header=F)[,1]),scientific=T),\"\n\",sep=\"\"))"
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      fgrep -v CTRL | cut -f1 )
  done | paste -s
done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed > \
${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal.bed 