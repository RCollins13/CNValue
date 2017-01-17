#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Performs a Binomial enrichment test for a class of genomic annotations in
# a set of case CNVs versus a set of control CNVs.

#Usage statement
usage(){
cat <<EOF
usage: CNV_set_annoClass_bulk_test.sh [-h] [-m ADJUSTMENT] [-a ALTERNATE] [-z]
                                      CONTROLS CASES ELEMENTS OUTFILE

Binomial enrichment test of genomic annotation class versus case/control CNVs

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  ELEMENTS   path to input list of BED-style annotations to test. Input list must
             have two columns: annotation ID and full path to annotation BED file.
             Each annotation BED file must have at least three columns: chr, start 
             coordinate, and end coordinate.
  OUTFILE    output file

Optional arguments:
  -h  HELP          Show this help message and exit
  -m  ADJUSTMENT    Case CNV count will be divided by this value to account for 
                    differences in case/control distributions (default: ratio
                    of case and control median sizes)
  -a  ALTERNATE     Alternative hypothesis to be tested (default: greater; 
                    options: greater, less)
  -z  GZIP          Gzip output file
EOF
}

#Parse arguments
ALTERNATE="greater"
GZ=0
ADJUSTMENT="size"
while getopts ":m:a:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    m)
      ADJUSTMENT=${OPTARG}
      ;;
    a)
      ALTERNATE=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
ELEMENTS=$3
OUTFILE=$4

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${ELEMENTS} ]; then
  usage
  exit 0
fi

#Makes master tmpdir for all working files
TMPDIR=`mktemp -d`

#Unzip input control CNVs file if gzipped
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

#Get counts and median size of CNVs per group
CTRL_CNVs=$( cat ${CTRL} | wc -l )
CASE_CNVs=$( cat ${CASE} | wc -l )
if [ ${ADJUSTMENT} == "size" ]; then
  CTRL_SIZE=$( awk '{ print $3-$2 }' ${CTRL} | sort -k1,1n | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  CASE_SIZE=$( awk '{ print $3-$2 }' ${CASE} | sort -k1,1n | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  SIZE_MOD=$( Rscript -e "cat(${CASE_SIZE}/${CTRL_SIZE})" | fgrep -v WARNING )
else
  SIZE_MOD=${ADJUSTMENT}
fi

#Prepares output file header
echo -e "#Annotation\tControl_Hits\tControl_Misses\tCase_Hits\tCase_Misses\tFoldChange\tFoldChange_CImin\tFoldChange_CImax\tp\t\
Case_Hits_Adj\tCase_Misses_Adj\tFoldChange_Adj\tFoldChange_CImin_Adj\tFoldChange_CImax_Adj\tp_Adj" > ${OUTFILE}

#Loops over annotations to test
while read bID BINS; do
  #Unzip input bins file if gzipped and truncate to 4 columns (at most)
  GZI_BINS=0
  BIN=`mktemp`
  if [ $( file ${BINS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
    GZI_BINS=1
    zcat ${BINS} | cut -f1-4 > ${BIN}
  else
    cut -f1-4 ${BINS} > ${BIN}
  fi

  #Count number of CNVs, those that hit at least one element, and those that miss all elements
  CTRL_HITS=$( bedtools intersect -u -a ${CTRL} -b ${BIN} | wc -l )
  CTRL_MISS=$( bedtools intersect -v -a ${CTRL} -b ${BIN} | wc -l )
  CASE_HITS=$( bedtools intersect -u -a ${CASE} -b ${BIN} | wc -l )
  CASE_MISS=$( bedtools intersect -v -a ${CASE} -b ${BIN} | wc -l )

  #Calculates statistics
  p=$( Rscript -e "cat(binom.test(x=${CASE_HITS},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs},alternative=\"${ALTERNATE}\")\$p.value)" | fgrep -v WARNING )
  BernoulliObs=$( Rscript -e "cat(binom.test(x=${CASE_HITS},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$estimate)" | fgrep -v WARNING )
  BernoulliExp=$( Rscript -e "cat(binom.test(x=${CASE_HITS},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$null.value)" | fgrep -v WARNING )
  ObsExp=$( Rscript -e "cat(${BernoulliObs}/${BernoulliExp})" | fgrep -v WARNING )
  CImin=$( Rscript -e "cat(binom.test(x=${CASE_HITS},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$conf.int[1]/${BernoulliExp})" | fgrep -v WARNING )
  CImax=$( Rscript -e "cat(binom.test(x=${CASE_HITS},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$conf.int[2]/${BernoulliExp})" | fgrep -v WARNING )

  #Adjusts case hits and misses according to differences in median size
  CASE_HITS_ADJ=$( Rscript -e "cat(round(${CASE_HITS}/${SIZE_MOD},0))" | fgrep -v WARNING )
  CASE_MISS_ADJ=$(( ${CASE_CNVs} - ${CASE_HITS_ADJ} ))

  #Calculates size-adjusted statistics
  p_ADJ=$( Rscript -e "cat(binom.test(x=${CASE_HITS_ADJ},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs},alternative=\"${ALTERNATE}\")\$p.value)" | fgrep -v WARNING )
  BernoulliObs_ADJ=$( Rscript -e "cat(binom.test(x=${CASE_HITS_ADJ},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$estimate)" | fgrep -v WARNING )
  BernoulliExp_ADJ=$( Rscript -e "cat(binom.test(x=${CASE_HITS_ADJ},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$null.value)" | fgrep -v WARNING )
  ObsExp_ADJ=$( Rscript -e "cat(${BernoulliObs_ADJ}/${BernoulliExp_ADJ})" | fgrep -v WARNING )
  CImin_ADJ=$( Rscript -e "cat(binom.test(x=${CASE_HITS_ADJ},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$conf.int[1]/${BernoulliExp_ADJ})" | fgrep -v WARNING )
  CImax_ADJ=$( Rscript -e "cat(binom.test(x=${CASE_HITS_ADJ},n=${CASE_CNVs},p=${CTRL_HITS}/${CTRL_CNVs})\$conf.int[2]/${BernoulliExp_ADJ})" | fgrep -v WARNING )

  #Prints result
  echo -e "${bID}\t${CTRL_HITS}\t${CTRL_MISS}\t${CASE_HITS}\t${CASE_MISS}\t${ObsExp}\t${CImin}\t${CImax}\t${p}\t\
${CASE_HITS_ADJ}\t${CASE_MISS_ADJ}\t${ObsExp_ADJ}\t${CImin_ADJ}\t${CImax_ADJ}\t${p_ADJ}"
done < ${ELEMENTS} >> ${OUTFILE}

#Gzip output if optioned
if [ ${GZ} -gt 0 ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
rm -rf ${TMPDIR}
