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

#Usage statement
usage(){
cat <<EOF
usage: geneSet_burdenTest_batch.sh [-h] [-N TIMES] [-U UNIVERSE] [-W WHOLE GENE]
                                   [-A ALLOSOMES] [-H OVERRIDE] [-o OUTDIR] [-f]
                                   CONTROLS CASES GENESET GTF

Script to test CNV burden at gene sets in batch mode

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  GENESET    path to gene set input file to be tested. Must have one column 
             with gene symbols
  GTF        path to GTF input file. Exons or gene boundaries will be extracted
             from this file

Optional arguments:
  -h  HELP          Show this help message and exit
  -U  UNIVERSE      List of gene symbols to consider as the background 
                    geneset (default: all genes)
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
OUTFILE=/dev/stdout
UNIVERSE=ALL
WG=0
ALLO=0
OVER=0
PREFIX=""
FORCE=0
while getopts ":N:U:W:A:H:p:o:hf" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    U)
      UNIVERSE=${OPTARG}
      ;;
    W)
      WG=1
      ;;
    A)
      ALLO=1
      ;;
    H)
      OVER=1
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
GENESET=$3
GTF=$4

###MUST ADD: #GET PATH TO RCNVMAP BIN SUBDIRECTORY
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${GENESET} ] || [ -z ${GTF} ]; then
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

# #Iterate over list of annotations and run burden tests
# while read NAME ANNO; do
#   if [ -e ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt ]; then
#     if [ ${FORCE} -eq 0 ]; then
#       echo "OUTPUT FILE FOR ${NAME} FOUND; SKIPPING"
#     else
#       echo "STARTING ${NAME}"
#       ${BIN}/geneSet_permutation_test.sh -q -N ${TIMES} -U ${UNIVERSE} -L ${NAME} \
#       -o ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt \
#       ${CONTROLS} ${CASES} ${ANNO} ${GENOME}
#     fi
#   else
#     echo "STARTING ${NAME}"
#     ${BIN}/annoSet_permutation_test.sh -q -N ${TIMES} -x ${EXCLUDE} -L ${NAME} \
#     -o ${OUTDIR}/${PREFIX}.${NAME}.CNV_burden_results.txt \
#     ${CONTROLS} ${CASES} ${ANNO} ${GENOME}
#   fi
# done < ${LIST}
