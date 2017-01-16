#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Shuffles CNV case/control labels

#Exit on error
set -e

#Usage statement
usage(){
cat <<EOF
usage: shuffle_CNV_labels.sh [-h] [-s] [-N times] [-o OUTDIR] [-p prefix] CONTROLS CASES TBRs

Shuffles CNV case/control labels

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  TBRs       path to TBR input file. Must have at least four columns: chr, 
             TBR min, TBR max, and TBR ID

Optional arguments:
  -h  HELP          Show this help message and exit
  -s  SAVE CNVs     Save shuffled CNV files (default: do not save)
  -N  TIMES         Number of permutations to perform (default: 100)
  -o  OUTDIR        Output directory (default: pwd)
  -p  PREFIX        Prefix for all output files (default: TBRden_shuffle)
EOF
}

#Parse arguments
OUTDIR=`pwd`
PREFIX="TBRden_shuffle"
SAVE=0
N=100
while getopts ":o:p:N:sh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTDIR=${OPTARG}
      ;;
    p)
      PREFIX=${OPTARG}
      ;;
    s)
      SAVE=1
      ;;
    N)
      N=${OPTARG}
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
TBRs=$3

#Set TBRden directory path
TBRden_bin="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ]; then
  usage
  exit 0
fi

#Checks to make sure OUTDIR exists; if not, attempts to create
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Unzip input control CNV file if gzipped
GZI_CTRL=0
if [ $( file ${CONTROLS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_CTRL=1
  CTRL=`mktemp`; mv ${CTRL} ${CTRL}.gz; CTRL=${CTRL}.gz
  cp ${CONTROLS} ${CTRL}
  gunzip ${CTRL}
  CTRL=$( echo "${CTRL}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  CTRL=CONTROLS
fi

#Unzip input case CNV file if gzipped
GZI_CASE=0
if [ $( file ${CASES} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_CASE=1
  CASE=`mktemp`; mv ${CASE} ${CASE}.gz; CASE=${CASE}.gz
  cp ${CASES} ${CASE}
  gunzip ${CASE}
  CASE=$( echo "${CASE}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  CASE=CASES
fi

#Count number of case & control CNVs
NCASE=$( fgrep -v "#" ${CASE} | wc -l )
NCTRL=$( fgrep -v "#" ${CTRL} | wc -l )

#Repeat for number of permutations specified by user
for i in $( seq 1 ${N} ); do

  #Create output directory for shuffle i
  mkdir ${OUTDIR}/${PREFIX}_${i}

  #Pool case & control CNVs and shuffle randomly
  ALL=`mktemp`
  cat <( fgrep -v "#" ${CASE} ) <( fgrep -v "#" ${CTRL} ) | cut -f1-4 | shuf > ${ALL}

  #Subset shuffled cases & controls
  SCASE=${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CASE_CNVs.bed
  SCTRL=${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CONTROL_CNVs.bed
  head -n ${NCASE} ${ALL} | sort -Vk1,1 -k2,2n -k3,3n > ${SCASE}
  tail -n ${NCTRL} ${ALL} | sort -Vk1,1 -k2,2n -k3,3n > ${SCTRL}
  gzip ${SCASE} ${SCTRL}
  SCASE=${SCASE}.gz
  SCTRL=${SCTRL}.gz

  #Run CNV pileups on shuffled CNVs
  ${TBRden_bin}/TBRden_pileup.sh -z -o ${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CASE_CNVs.pileup.bed.gz ${SCASE} ${TBRs}
  ${TBRden_bin}/TBRden_pileup.sh -z -o ${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CONTROL_CNVs.pileup.bed.gz ${SCTRL} ${TBRs}

  #Calculate burden from shuffled CNVs
  ${TBRden_bin}/TBRden_test.R \
  ${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CONTROL_CNVs.pileup.bed.gz \
  ${OUTDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CASE_CNVs.pileup.bed.gz \
  ${OUTDIR}/${PREFIX}_${i} ${PREFIX}_${i}

  #Clean up
  if [ SAVE -eq 0 ]; then
    rm ${SCASE} ${SCTRL}
  fi
done
