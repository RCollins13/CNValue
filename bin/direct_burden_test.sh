#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Direct test of CNV burden at a given loci by localized CNV shifting

#Usage statement
usage(){
cat <<EOF
usage: direct_burden_test.sh [-h] [-d DIST] [-N TIMES] [-e EXONS] [-c/-n] [-t TAIL]
                             [-p prefix] [-z] CONTROLS CASES ELEMENTS OUTDIR

Direct test of CNV burden at a given loci by localized CNV shifting

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  ELEMENTS   path to input file of BED-style bins to test. Must have at least 
             three columns: chr, start coordinate, and end coordinate. May 
             optionally have a fourth column for element ID.
  OUTDIR     output directory

Optional arguments:
  -h  HELP          Show this help message and exit
  -d  DIST          Maximum distance to shift each CNV, as fraction of 
                    CNV size (default: 5.0)
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -e  EXONS         File of exons to use for coding/noncoding filter
  -c  CODING        Flag to only consider coding (i.e. exon-overlapping) CNVs
  -n  NONCODING     Flag to only consider noncoding (i.e. not exon-overlapping) CNVs
  -t  TAIL          Assess "upper" or "lower" tail for p-value reporting
  -p  PREFIX        Prefix for all output files (default: TBRden_shuffle)
  -z  GZIP          Gzip output
EOF
}

#Parse arguments
DIST=5
TIMES=1000
CODING=0
NONCODING=0
EXONS=0
TAIL="upper"
PREFIX="TBRden_direct_burden_test"
GZ=0
while getopts ":d:N:e:t:p:zcnh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    d)
      DIST=${OPTARG}
      ;;
    b)
      BUFFER=${OPTARG}
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    e)
      EXONS=${OPTARG}
      ;;
    t)
      TAIL=${OPTARG}
      ;;
    c)
      CODING=1
      ;;
    n)
      NONCODING=1
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
CONTROLS=$1
CASES=$2
BINS=$3
OUTDIR=$4

#Check that both coding and noncoding flags aren't set
if [ ${CODING} -eq 1 ] && [ ${NONCODING} -eq 1 ]; then
  echo -e "\nERROR: SPECIFY EITHER -c OR -n, BUT NOT BOTH\n"
  usage
  exit 0
fi

#Check that exon file exists if either coding or noncoding flag is set
if [ ${CODING} -eq 1 ] || [ ${NONCODING} -eq 1 ]; then
  if ! [ -e ${EXONS} ]; then
    echo -e "\nERROR: INVALID EXON FILE (REQUIRED FOR EITHER -c OR -n)\n"
    usage
    exit 0
  fi
fi

#Check that output directory exists; attempt to create if not
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Set TBRden directory path
TBRden_bin="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${BINS} ]; then
  usage
  exit 0
fi

#Makes master tmpdir for all working files
TMPDIR=`mktemp -d`

#Unzip input bins file if gzipped
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

#Unzip input bins file if gzipped
GZI_BINS=0
if [ $( file ${BINS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_BINS=1
  BIN=`mktemp`; mv ${BIN} ${BIN}.gz; BIN=${BIN}.gz
  cp ${BIN} ${BIN}
  gunzip ${BIN}
  BIN=$( echo "${CASE}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  BIN=${BINS}
fi

#Calculates observed dCNVs for all sites
OBSERVED=`mktemp`
paste <( bedtools intersect -c -a ${BIN} -b ${CASE} | awk '{ print $NF }' ) \
<( bedtools intersect -c -a ${BIN} -b ${CTRL} | awk '{ print $NF }' ) | \
awk -v OFS="\t" '{ print $1-$2 }' | paste <( fgrep -v "#" ${BINS} ) - > ${OBSERVED}

#Get counts of number of case & control CNVs
nCASE=$( fgrep -v "#" ${CASE} | wc -l )
nCTRL=$( fgrep -v "#" ${CTRL} | wc -l )

#Make temporary files for iterating permutations
DIRECTION=`mktemp`
CASE_SHUF=`mktemp`
CTRL_SHUF=`mktemp`

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  #Print status
  echo "Beginning permutation ${i} of ${TIMES}"
  
  #Create output directory for shuffle i
  mkdir ${TMPDIR}/${PREFIX}_${i}

  #Shift case CNVs by half their size
  Rscript -e "write.table(as.data.frame(sample(seq(-${DIST},${DIST},0.01),${nCASE},replace=T)),\
              \"${DIRECTION}\",row.names=F,col.names=F,quote=F)"
  if [ ${CODING} -eq 1 ]; then
    paste <( fgrep -v "#" ${CASE} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' | \
    bedtools intersect -wa -u -a - -b ${EXONS} > ${CASE_SHUF}
  elif [ ${NONCODING} -eq 1 ]; then
    paste <( fgrep -v "#" ${CASE} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' | \
    bedtools intersect -wa -v -a - -b ${EXONS} > ${CASE_SHUF}
  else
    paste <( fgrep -v "#" ${CASE} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' > ${CASE_SHUF}
  fi

  #Shift control CNVs by half their size
  Rscript -e "write.table(as.data.frame(sample(seq(-${DIST},${DIST},0.01),${nCTRL},replace=T)),\
              \"${DIRECTION}\",row.names=F,col.names=F,quote=F)"
  if [ ${CODING} -eq 1 ]; then
    paste <( fgrep -v "#" ${CTRL} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' | \
    bedtools intersect -wa -u -a - -b ${EXONS} > ${CTRL_SHUF}
  elif [ ${NONCODING} -eq 1 ]; then
    paste <( fgrep -v "#" ${CTRL} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' | \
    bedtools intersect -wa -v -a - -b ${EXONS} > ${CTRL_SHUF}
  else
    paste <( fgrep -v "#" ${CTRL} ) ${DIRECTION} | awk -v OFS="\t" \
    '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
    awk -v OFS="\t" '{ if ($2>=0) print }' > ${CTRL_SHUF}
  fi

  #Run CNV pileups on shuffled CNVs
  paste <( bedtools intersect -c -a ${BIN} -b ${CASE_SHUF} | awk '{ print $NF }' ) \
  <( bedtools intersect -c -a ${BIN} -b ${CTRL_SHUF} | awk '{ print $NF }' ) | \
  awk -v OFS="\t" '{ print $1-$2 }' | paste <( fgrep -v "#" ${BINS} ) - > \
  ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}_pileups.bed
done

#Collect results
dCNV_RES=`mktemp`
for i in $( seq 1 ${TIMES} ); do
  awk '{ print $NF }' ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}_pileups.bed | paste -s
done > ${dCNV_RES}

#Run statistical modeling & print results to outfile
${TBRden_bin}/TBRden_direct_burden_test.R \
${OBSERVED} \
${dCNV_RES} \
${TAIL} \
brown \
${OUTDIR} \
${PREFIX}

#Clean up
rm -rf ${TMPDIR}
