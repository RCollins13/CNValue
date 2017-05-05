#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to curate gene lists for each GO term

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read argument
list=$1

#####Run
while read GO term; do
  echo ${term}
  nocolon=$( echo "${GO}" | sed 's/\:/_/g' )
  fgrep -w ${GO} ${WRKDIR}/data/misc/GO/goa_human.gaf | cut -f3 | sort | uniq > \
  ${WRKDIR}/data/misc/GO/all_GO_gene_lists/${nocolon}_${term}.genes.list
  if ! [ -s ${WRKDIR}/data/misc/GO/all_GO_gene_lists/${nocolon}_${term}.genes.list ]; then
    rm ${WRKDIR}/data/misc/GO/all_GO_gene_lists/${nocolon}_${term}.genes.list
  fi
done < ${list}