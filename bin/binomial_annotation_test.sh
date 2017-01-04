#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Binomial test of genomic annotation enrichment against a set of loci of interest

#Usage statement
usage(){
cat <<EOF
usage: binomial_annotation_test.sh [-h] [-N TIMES] [-x EXCLUDE] [-a ALTERNATIVE]
                                   [-p prefix] [-z] LOCI ALL_BINS ANNO OUTDIR 

Binomial test of genomic annotation enrichment against a set of loci of interest

Positional arguments:
  LOCI       path to BED file of loci to test
  ALL_BINS   path to BED file of all bins eligible to be selected during permutation
  ANNO       path to BED file of annotations to consider
  OUTDIR     output directory

Optional arguments:
  -h  HELP          Show this help message and exit
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -x  EXCLUDE       BED intervals used to blacklist annotation file
  -a  ALTERNATIVE   Alternative hypothesis to test; either "greater" or "less"
  -p  PREFIX        Prefix for all output files (default: TBRden_binomial_annotation_test)
  -z  GZIP          Gzip output
EOF
}

#Parse arguments
TIMES=1000
EXCLUDE=0
ALT="upper"
PREFIX="TBRden_binomial_annotation_test"
GZ=0
while getopts ":N:x:a:p:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    x)
      EXCLUDE=${OPTARG}
      ;;
    a)
      ALT=${OPTARG}
      ;;
    p)
      PREFIX=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
LOCI=$1
ALL_BINS=$2
ANNO=$3
OUTDIR=$4

#Check that exclusion file exists if specified
if [ ${EXCLUDE} != "0" ] && ! [ -e ${EXCLUDE} ]; then
  echo -e "\nERROR: EXCLUSION FILE NOT FOUND\n"
  usage
  exit 0
fi

#Check that output directory exists; attempt to create if not
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Set TBRden directory path
TBRden_bin="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Check for required input
if [ -z ${LOCI} ] || [ -z ${ALL_BINS} ] || [ -z ${ANNO} ]; then
  usage
  exit 0
fi

#Makes master tmpdir for all working files
TMPDIR=`mktemp -d`

#Unzip input loci file if gzipped
GZI_LOCI=0
if [ $( file ${LOCI} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_LOCI=1
  LOCI1=`mktemp`; mv ${LOCI1} ${LOCI1}.gz; LOCI1=${LOCI1}.gz
  cp ${LOCI} ${LOCI1}
  gunzip ${LOCI1}
  LOCI=$( echo "${LOCI1}" | sed 's/\.gz/\t/g' | cut -f1 )
fi

#Trim loci file to four columns if larger
TEST_LOCI=`mktemp`
cut -f1-4 ${LOCI} > ${TEST_LOCI}

#Unzip input all bins file if gzipped
GZI_BINS=0
if [ $( file ${ALL_BINS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_BINS=1
  BINS=`mktemp`; mv ${BINS} ${BINS}.gz; BINS=${BINS}.gz
  cp ${ALL_BINS} ${BINS}
  gunzip ${BINS}
  BINS=$( echo "${BINS}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  BINS=${ALL_BINS}
fi

#Trim bins file to four columns if larger
TEST_BINS=`mktemp`
cut -f1-4 ${BINS} > ${TEST_BINS}

#Blacklists annotation file if optioned
CLEAN_ANNO=`mktemp`
if [ ${EXCLUDE} == "0" ]; then
  cp ${ANNO} ${CLEAN_ANNO}
else
  bedtools intersect -wa -a ${ANNO} -b ${EXCLUDE} > ${CLEAN_ANNO}
fi

#Calculates observed annotation overlap for all sites
OBSERVED=`mktemp`
bedtools intersect -c -a ${TEST_LOCI} -b ${CLEAN_ANNO} > ${OBSERVED}

#Gets number and sizes of test loci
SIZES=`mktemp`
fgrep -v "#" ${TEST_LOCI} | awk '{ print $3-$2 }' > ${SIZES}
nLOCI=$( fgrep -v "#" ${TEST_LOCI} | wc -l )

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  #Generate new bins & intersect with annotations
  hits=$( fgrep -v "#" ${TEST_BINS} | shuf | head -n${nLOCI} | paste - ${SIZES} | \
  awk -v OFS="\t" '{ print $1, $2, $2+$NF }' | bedtools intersect -c -a - -b ${CLEAN_ANNO} | \
  awk '{ if ($NF>0) print $0 }' | wc -l )

  #Print results
  echo "${hits}"
done > ${TMPDIR}/perm_results.txt

#Run helper R script to calculate P value and plot results

#Run statistical modeling & print results to outfile
${TBRden_bin}/TBRden_run_binomialHelper.R \
$( awk '{ if ($NF>0) print $0 }' ${OBSERVED} | wc -l ) \
${nLOCI} \
${TMPDIR}/perm_results.txt \
${ALT} \
${OUTDIR} \
${PREFIX}

#Clean up
rm -rf ${TMPDIR}
