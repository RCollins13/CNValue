#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Script to iterate over a list of gene sets and test each for
# CNV burden from a specified set of CNVs in cases and controls

# #Testing dev parameters
# TIMES=100
# UNIVERSE=${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list
# WG=0
# ALLO=0
# OVER=${WRKDIR}/data/misc/exons_boundaries_dictionary/
# PREFIX="geneSetTest"
# OUTDIR=${TMPDIR}/batchGeneSetTests/
# FORCE=1
# CONTROLS=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.E3.GRCh37.all.bed.gz
# CASES=${WRKDIR}/data/CNV/CNV_MASTER/NDD/NDD.DEL.E3.GRCh37.all.bed.gz
# LIST=/data/talkowski/Samples/rCNVmap/bin/rCNVmap/misc/master_gene_sets.list
# GTF=${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf
# BIN=/data/talkowski/Samples/rCNVmap/bin/rCNVmap/bin/

#Usage statement
usage(){
cat <<EOF
usage: geneSet_burdenTest_batch.sh [-h] [-N TIMES] [-W WHOLE GENE] [-A ALLOSOMES]
                                   [-H OVERRIDE] [-p PREFIX] [-o OUTDIR] [-f]
                                   CONTROLS CASES LIST GTF

Script to test CNV burden at gene sets in batch mode

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  LIST       three-column, tab-delimmed list of gene sets to test.
             First column: gene set name; second column: full path to gene set;
             Third column: full path to universe gene set
  GTF        path to GTF input file. Exons or gene boundaries will be extracted
             from this file. Can be overridden with -H

Optional arguments:
  -h  HELP          Show this help message and exit
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -W  WHOLE GENE    Restrict analysis to CNVs that span the entire gene
                    (default: count any exonic overlap)
  -A  ALLOSOMES     Include allosomes in analyses (default: false)
  -H  OVERRIDE      Path to correctly formatted exons & boundaries directory
                    Note: do not use this option unless you know what you're doing
  -p  PREFIX        String to add to the front of all output files
  -o  OUTDIR        Output directory (default: current directory)
  -f  FORCE         Overwrite output even if file with same name exists
EOF
}

#Parse arguments
TIMES=1000
WG=0
ALLO=0
OVER=0
PREFIX="geneSetTest"
OUTDIR=`pwd`
FORCE=0
while getopts ":N:W:A:H:p:o:hf" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    W)
      WG=1
      ;;
    A)
      ALLO=1
      ;;
    H)
      OVER=${OPTARG}
      ;;
    p)
      PREFIX=${OPTARG}
      ;;
    o)
      OUTDIR=${OPTARG}
      ;;
    f)
      FORCE=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
LIST=$3
GTF=$4

###MUST ADD: #GET PATH TO RCNVMAP BIN SUBDIRECTORY
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${LIST} ] || [ -z ${GTF} ]; then
  usage
  exit 0
fi

#Attempts to create $OUTDIR if it doesn't exist
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Exit if $OUTDIR doesn't exist
if ! [ -e ${OUTDIR} ]; then
  usage
  echo -e "\nERROR: OUTDIR DOES NOT EXIST"
  exit 0
fi

#Write list of based options
opts=$( echo -e "-q -N ${TIMES}" )
if [ ${WG} -gt 0 ]; then
  opts="${opts} -W"
fi
if [ ${ALLO} -gt 0 ]; then
  opts="${opts} -A"
fi
if [ ${OVER}!="0" ]; then
  opts="${opts} -H ${OVER}"
fi

#Iterate over list of annotations and run burden tests
while read NAME GENESET UNIVERSE; do
  if [ -e ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt ]; then
    if [ ${FORCE} -eq 0 ]; then
      echo "OUTPUT FILE FOR ${NAME} FOUND; SKIPPING"
    else
      echo "STARTING ${NAME}"
      ${BIN}/geneSet_permutation_test.sh ${opts} -L ${NAME} -U ${UNIVERSE} \
      -o ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt \
      ${CONTROLS} ${CASES} ${GENESET} ${GTF}
    fi
  else
    echo "STARTING ${NAME}"
    ${BIN}/geneSet_permutation_test.sh ${opts} -L ${NAME} -U ${UNIVERSE} \
    -o ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt \
    ${CONTROLS} ${CASES} ${GENESET} ${GTF}
  fi
done < ${LIST}
