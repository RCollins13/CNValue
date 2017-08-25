#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects data required for modeling per-element rCNV burden scores

#Testing dev parameters
# export WRKDIR=/data/talkowski/Samples/rCNVmap
# source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh
# PREFIX="BRAIN_MASTER.UCNE"
# OUTFILE=${TMPDIR}/annoset_test.out
# WHOLE=0
# ALLO=0
# GZIP=1
# QUIET=0
# CONTROLS=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.E4.GRCh37.haplosufficient.bed.gz
# CASES=${WRKDIR}/data/CNV/CNV_MASTER/NDD/NDD.DEL.E4.GRCh37.haplosufficient.bed.gz
# BED=${WRKDIR}/data/master_annotations/noncoding/UCNEs.elements.bed
# REF=${h37}
# BIN=${WRKDIR}/bin/rCNVmap/bin/

#Usage statement
usage(){
cat <<EOF
usage: gather_annoScore_data.sh [-h] [-p PREFIX] [-W WHOLE] [-A ALLOSOMES] [-z]
                                [-o OUTFILE] [-q QUIET] CONTROLS CASES GTF REF

Collects data required for modeling per-element rCNV burden scores

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  BED        path to BED input file of elements to test
  REF        path to reference .fa 

Optional arguments:
  -h  HELP          Show this help message and exit
  -p  PREFIX        Prefix to name each element in test set
  -W  WHOLE ELEMENT Restrict analysis to CNVs that span the entire element
                    (default: count any overlap)
  -A  ALLOSOMES     Include allosomes in analyses (default: false)
  -z  GZIP          Gzip output file (default: false)
  -o  OUTFILE       Output file (default: /dev/stdout)
  -q  QUIET         Suppresses (some) standard output
EOF
}

#Parse arguments
PREFIX="ELEMENT"
OUTFILE=/dev/stdout
WHOLE=0
ALLO=0
GZIP=0
QUIET=0
while getopts ":p:WAzo:qh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    W)
      WHOLE=1
      ;;
    p)
      PREFIX=1
      ;;
    A)
      ALLO=1
      ;;
    z)
      GZIP=1
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    q)
      QUIET=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
BED=$3
REF=$4

#Get path to rCNVmap bin
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${BED} ] || [ -z ${REF} ]; then
  usage
  exit 0
fi

#Makes master tmpdir for all working files
TMPDIR=`mktemp -d`

#Unzip input control CNV file if gzipped
GZI_CTRL=0
if [ $( file ${CONTROLS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_CTRL=1
  CTRL=`mktemp`; mv ${CTRL} ${CTRL}.gz; CTRL=${CTRL}.gz
  cp ${CONTROLS} ${CTRL}
  gunzip ${CTRL}
  CTRL=$( echo "${CTRL}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  CTRL=${CONTROLS}
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
  CASE=${CASES}
fi

#Restricts CNVs to autosomes unless optioned
if [ ${ALLO} -eq 0 ]; then
  grep -e '^[0-9]\|^chr[0-9]' ${CTRL} > ${CTRL}2; mv ${CTRL}2 ${CTRL}
  grep -e '^[0-9]\|^chr[0-9]' ${CASE} > ${CASE}2; mv ${CASE}2 ${CASE}
fi

#Format elements to test
ELEMENTS=`mktemp`
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::HARMONIZING INPUT DATASETS..."
fi
if [ ${ALLO} -eq 0 ]; then
  grep -e '^[0-9]\|^chr[0-9]' ${BED} | awk -v PREFIX=${PREFIX} -v OFS="\t" \
  '{ print $1, $2, $3, PREFIX"_"NR }' > ${ELEMENTS}
else
  awk -v PREFIX=${PREFIX} -v OFS="\t" \
  '{ print $1, $2, $3, PREFIX"_"NR }' ${BED} > ${ELEMENTS}
fi

#Annotates GC content of elements
GC_PROFILES=`mktemp`
bedtools nuc -fi ${REF} -bed ${ELEMENTS} | fgrep -v "#" | cut -f4,6 | \
sort -k1,1 > ${GC_PROFILES}

#Isolate unique CNV-element pairs
CNV_ELEMENT_PAIRS_CASE=`mktemp`
CNV_ELEMENT_PAIRS_CTRL=`mktemp`
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::OVERLAPPING CNVS AND ELEMENTS..."
fi
if [ ${WHOLE} -eq 1 ]; then
  #Case, whole-anno CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${ELEMENTS} -b ${CASE} | \
  awk -v OFS="\t" '{ print $8, $4 }' | \
  sort -k1,1 -k2,2 | uniq > ${CNV_ELEMENT_PAIRS_CASE}
  #Control, whole-anno CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${ELEMENTS} -b ${CTRL} | \
  awk -v OFS="\t" '{ print $8, $4 }' | \
  sort -k1,1 -k2,2 | uniq > ${CNV_ELEMENT_PAIRS_CTRL}
else
  #Case, exonic CNVs
  bedtools intersect -wa -wb -a ${CASE} -b ${ELEMENTS} | \
  awk -v OFS="\t" '{ print $4, $NF }' | sort -k1,1 -k2,2 | \
  uniq > ${CNV_ELEMENT_PAIRS_CASE}
  #Control, exonic CNVs
  bedtools intersect -wa -wb -a ${CTRL} -b ${ELEMENTS} | \
  awk -v OFS="\t" '{ print $4, $NF }' | sort -k1,1 -k2,2 | \
  uniq > ${CNV_ELEMENT_PAIRS_CTRL}
fi

#Cut list of all element IDs
ELEMENT_IDs=`mktemp`
cut -f4 ${ELEMENTS} > ${ELEMENT_IDs}

#Calculate baseline dCNV & dElement info
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::RUNNING CNV WEIGHTING MODEL..."
fi
${BIN}/gather_annoScore_data.helper.R \
${CNV_ELEMENT_PAIRS_CASE} \
${CNV_ELEMENT_PAIRS_CTRL} \
${ELEMENT_IDs} \
${TMPDIR}/

#Write header to output file
echo -e "#chr\tstart\tend\telement_ID\telement_size\tGC\tcase_CNV\tcontrol_CNV\tcase_CNV_weighted\tcontrol_CNV_weighted" > ${OUTFILE}

#Collect final output
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::SUMMARIZING RESULTS PER ELEMENT..."
fi
nElements=$( cat ${ELEMENTS} | wc -l )
if [ $( cat ${GC_PROFILES} | wc -l ) -eq ${nElements} ] && \
   [ $( cat ${TMPDIR}/CASE.CNVsPerAnno.txt | wc -l ) -eq ${nElements} ] && \
   [ $( cat ${TMPDIR}/CTRL.CNVsPerAnno.txt | wc -l ) -eq ${nElements} ] && \
   [ $( cat ${TMPDIR}/CASE.weightedCNVsPerAnno.txt | wc -l ) -eq ${nElements} ] && \
   [ $( cat ${TMPDIR}/CTRL.weightedCNVsPerAnno.txt | wc -l ) -eq ${nElements} ]; then
  paste <( awk -v OFS="\t" '{ print $1, $2, $3, $4, $3-$2 }' ${ELEMENTS} ) \
    <( cut -f2 ${GC_PROFILES} ) \
    <( cut -f2 ${TMPDIR}/CASE.CNVsPerAnno.txt ) \
    <( cut -f2 ${TMPDIR}/CTRL.CNVsPerAnno.txt ) \
    <( cut -f2 ${TMPDIR}/CASE.weightedCNVsPerAnno.txt ) \
    <( cut -f2 ${TMPDIR}/CTRL.weightedCNVsPerAnno.txt ) >> \
    ${OUTFILE}
  else
    echo "ERROR: some elements missing from pre-output files! Exiting..."
    exit 0
  fi
fi

#Gzip, if optioned
if [ ${GZIP} -eq 1 ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
rm -rf ${TMPDIR}
