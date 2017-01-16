#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Runs intersection of a CNV dataset versus a TAD boundary file

#Usage statement
usage(){
cat <<EOF
usage: TBRden_pileup.sh [-h] [-z] [-o OUTFILE] CNVs TADs

Runs intersection of a CNV dataset versus a TAD boundary file

Positional arguments:
  CNVs     path to CNV input file. Must have at least three columns: chr, 
           start, end
  TADs     path to TAD input file. Must have at least seven columns: chr, 
           TBR1 min, TBR1 max, chr, TBR2 min, TBR2 max, TAD ID. Coordinates
           for TBRs insulating adjacent TADs must match exactly. No overlap 
           between non-identical TBRs is allowed. 

Optional arguments:
  -h  HELP          Show this help message and exit
  -z  GZIP          Gzip output file
  -o  OUTFILE       Output file (default: stdout)
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
GZ=0
while getopts ":o:m:zsh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
CNVs=$1
TADs=$2

#Check for required input
if [ -z ${CNVs} ] || [ -z ${TADs} ]; then
  usage
  exit 0
fi

#Check for gzip optioned only if output file specified
if [ ${OUTFILE} == "/dev/stdout" ] && [ ${GZ} == 1 ]; then
  echo -e "\nOUTFILE required for zip-compressed output"
  usage
  exit 0
fi

#Scrub ".gz" from output filename if provided by user
if [ ${GZ} == 1 ] && [ ${OUTFILE: -3} == ".gz" ]; then
  OUTFILE=$( echo "${OUTFILE}" | sed 's/\.gz//g' )
fi

#Unzip input TAD file if gzipped
GZI=0
if [ $( file ${TADs} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI=1
  TADI=`mktemp`; mv ${TADI} ${TADI}.gz; TADI=${TADI}.gz
  cp ${TADs} ${TADI}
  gunzip ${TADI}
  TADs=$( echo "${TADI}" | sed 's/\.gz/\t/g' | cut -f1 )
fi

#Link each left TBR to all overlapping right TBRs, and vice versa
#Merge identical TBRs used across multiple TADs
#Report flanking TAD coordinates to test for each TBR
ALLINTS=`mktemp`
fgrep -v "#" ${TADs} | awk -v OFS="\t" \
  '{ print $1, $2, $3, $7"_TBR_L\n"$1, $3, $5, $7"_M\n"$1, $5, $6, $7"_TBR_R" }' | \
  sort -Vk1,1 -k2,2n -k3,3n > ${ALLINTS}
TESTINTS=`mktemp`
while read chr start end TBRID; do
  side=$( echo ${TBRID} | sed 's/_/\t/g' | awk '{ print $NF }' )
  TBR_size=$((${end}-${start}))
  #Get coordinates of native flanking TAD
  native=$( echo ${TBRID} | sed 's/_TBR/\t/g' | awk '{ print $1"_M" }' )
  native_chr=$( fgrep -w ${native} ${ALLINTS} | cut -f1 )
  native_start=$( fgrep -w ${native} ${ALLINTS} | cut -f2 )
  native_end=$( fgrep -w ${native} ${ALLINTS} | cut -f3 )
  #Check for matching TBR
  if [ ${side} == "L" ]; then
    match=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( fgrep "_R" ${ALLINTS} ) | \
      awk '{ print $7 }' | paste -s -d, )
  else
    match=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( fgrep "_L" ${ALLINTS} ) | \
      awk '{ print $7 }' | paste -s -d, )
  fi
  #If matching TBR doesn't exist, report NAs for all fields
  if [ -z ${match} ]; then
    opposite=NA
    opposite_chr=NA
    opposite_start=NA
    opposite_end=NA
  else
    #Take smallest TAD of matched flanking TAD (after restricting TBR interval)
    #Smallest matching TAD must be at least as large as TBR
    opposite=$( echo "${match}" | sed -e 's/,/\n/g' -e 's/_TBR/\t/g' | \
      awk '{ print $1"_M" }' | fgrep -wf - ${ALLINTS} | \
      bedtools subtract -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | \
      awk -v OFS="\t" -v TBR_size=${TBR_size} '{ if ($3-$2>=TBR_size) print $0, $3-$2 }' | \
      sort -nk5,5 | cut -f4 | head -n1 )
    #If all matches are too small, report NAs for all fields
    if [ -z ${opposite} ]; then
      opposite=NA
      opposite_chr=NA
      opposite_start=NA
      opposite_end=NA
    else
      opposite_chr=$( fgrep -w ${opposite} ${ALLINTS} | cut -f1 )
      opposite_start=$( fgrep -w ${opposite} ${ALLINTS} | \
        bedtools subtract -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | cut -f2 )
      opposite_end=$( fgrep -w ${opposite} ${ALLINTS} | \
        bedtools subtract -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | cut -f3 )
    fi
  fi
  #Report TBR, left TAD, and right TAD
  if [ ${side} == "L" ]; then
    echo -e "${chr}\t${start}\t${end}\t${TBRID}\t\
    ${opposite_chr}\t${opposite_start}\t${opposite_end}\t${opposite}\t\
    ${native_chr}\t${native_start}\t${native_end}\t${native}"
  else
    echo -e "${chr}\t${start}\t${end}\t${TBRID}\t\
    ${native_chr}\t${native_start}\t${native_end}\t${native}\t\
    ${opposite_chr}\t${opposite_start}\t${opposite_end}\t${opposite}"
  fi
done < <( awk '$4 ~ /_L|_R/ { print $0 }' ${ALLINTS} ) > ${TESTINTS}

#Count number of CNVs over each interval being considered
#Note: score reports CNVs normalized by size (CNVs per kb); this value is *not* the
#same as the CNV density over the interval (i.e. average # of CNVs observed per base)
TMPCOV=`mktemp`
awk -v OFS="\t" '{ print $1, $2, $3, $4"\n"$5, $6, $7, $4"_"$8"\n"$9, $10, $11, $4"_"$12 }' ${TESTINTS} | \
awk '{ if ($1!="NA") print $0 }' | bedtools intersect -c -a - -b ${CNVs} | awk -v OFS="\t" \
  '{ print $1, $2, $3, $4, $5, 1000*$5/($3-$2) }' > ${TMPCOV}

#Consolidate output and report counts per TBR
echo -e "#chr\tstart\tend\tTBRs\tTBR_count\tTBR_score\tTAD_L\tTAD_L_count\tTAD_L_score\tTAD_R\tTAD_R_count\tTAD_R_score" > ${OUTFILE}
while read TBR_chr TBR_start TBR_end TBR_ID TADL_chr TADL_start TADL_end TADL_ID TADR_chr TADR_start TADR_end TADR_ID; do
  #Get count & score of TBR
  TBR_count=$( fgrep -w ${TBR_ID} ${TMPCOV} | cut -f5 )
  TBR_score=$( fgrep -w ${TBR_ID} ${TMPCOV} | cut -f6 )
  #Gather size and counts per flanking TAD
  if ! [ ${TADL_ID} == "NA" ]; then
    TADL_count=$( fgrep -w "${TBR_ID}_${TADL_ID}" ${TMPCOV} | cut -f5 )
    TADL_score=$( fgrep -w "${TBR_ID}_${TADL_ID}" ${TMPCOV} | cut -f6 )
  else
    TADL_count="NA"
    TADL_score="NA"
  fi
  if ! [ ${TADR_ID} == "NA" ]; then
    TADR_count=$( fgrep -w "${TBR_ID}_${TADR_ID}" ${TMPCOV} | cut -f5 )
    TADR_score=$( fgrep -w "${TBR_ID}_${TADR_ID}" ${TMPCOV} | cut -f6 )
  else
    TADR_count="NA"
    TADR_score="NA"
  fi
  #Print results
  echo -e "${TBR_chr}\t${TBR_start}\t${TBR_end}\t${TBR_ID}\t${TBR_count}\t${TBR_score}\t\
  ${TADL_ID}\t${TADL_count}\t${TADL_score}\t${TADR_ID}\t${TADR_count}\t${TADR_score}"
done < ${TESTINTS} >> ${OUTFILE}

#Gzip OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
if [ GZI == 1 ]; then
  rm -rf ${TMPI} ${TADs}
fi
rm -rf ${TMPCOV} ${ALLINTS} ${TESTINTS}
