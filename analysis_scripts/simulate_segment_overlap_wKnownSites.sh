#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Simulate expected overlap versus known sites

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh
VF=E2; filt=all

#####Read arguments
CNV=$1
ref=$2

#####Run
for i in $( seq 1 100 ); do
  bedtools shuffle -noOverlapping -seed ${i} \
  -excl ${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed \
  -i ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed \
  -g /data/talkowski/rlc47/src/GRCh37.genome | \
  bedtools coverage -b - -a ${WRKDIR}/data/misc/known_CNV_segments/${ref}.bed | \
  awk '{ if ($NF>=0.5) print $1 }' | wc -l
done > ${TMPDIR}/${CNV}.${ref}.simulated.txt
