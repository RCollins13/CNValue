#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to process & analyze all ABC EP connections from Lander lab

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Prep directories
for dir in ABC_EP ABC_EP/RLC_processed; do
  if ! [ -e ${WRKDIR}/data/misc/${dir} ]; then
    mkdir ${WRKDIR}/data/misc/${dir}
  fi
done

#####Rsync data
bsub -q filemove -sla miket_sc -J ABC_EP_rsync \
"rsync -ravz \
rcollins@gold.broadinstitute.org:/seq/lincRNA/RAP/ABC/171116_pred/* \
${WRKDIR}/data/misc/ABC_EP/"
