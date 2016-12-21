#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Runs intersection of a CNV dataset versus a TBR bed file

#Usage statement
usage(){
cat <<EOF
usage: TBRden_pileup.sh [-h] [-z] [-D distance] [-o OUTFILE] CNVs TBRs

Runs intersection of a CNV dataset versus a TAD boundary region (TBR) file

Positional arguments:
  CNVs     path to CNV input file. Must have at least three columns: chr, 
           start, end
  TBRs     path to TBR input file. Must have at least four columns: chr, 
           TBR min, TBR max, and TBR ID.

Optional arguments:
  -h  HELP          Show this help message and exit
  -D  DISTANCE      Flank distance for normalization, in bp (default 500,000)
  -z  GZIP          Gzip output file
  -o  OUTFILE       Output file (default: stdout)
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
GZ=0
DIST=500000
while getopts ":o:m:D:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    D)
      DIST=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
CNVs=$1
TBRs=$2

#Check for required input
if [ -z ${CNVs} ] || [ -z ${TBRs} ]; then
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
if [ $( file ${TBRs} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI=1
  TADI=`mktemp`; mv ${TADI} ${TADI}.gz; TADI=${TADI}.gz
  cp ${TBRs} ${TADI}
  gunzip ${TADI}
  TBRs=$( echo "${TADI}" | sed 's/\.gz/\t/g' | cut -f1 )
fi

#Make master list of all intervals needed to do pileup calculations
INTS=`mktemp`
fgrep -v "#" ${TBRs} | awk -v OFS="\t" -v D=${DIST} \
'{ print $1, $2-D, $2, $4"_A\n"$1, $2, $3, $4"_B\n"$1, $3, $3+D, $4"_C" }' | \
awk -v OFS="\t" '{ if ($2<=1) $2=1; print }' | awk -v OFS="\t" '{ if ($3<=2) $3=2; print }' > ${INTS}
# #Used for region-specific counts (not used in current version)
# fgrep -v "#" ${TBRs} | awk -v OFS="\t" -v D=${DIST} \
# '{ print $1, $2-D, $2, $4"_A\n"$1, $2, $3, $4"_B\n"$1, $3, $3+D, $4"_C\n"\
# $1, $2-D, $3, $4"_AB\n"$1, $2, $3+D, $4"_BC\n"$1, $2-D, $3+D, $4"_ABC" }' > ${INTS}

#Run overlap of master interval file vs CNVs
COUNTS=`mktemp`
bedtools intersect -c -a ${INTS} -b ${CNVs} > ${COUNTS}

#Write results to OUTFILE with header
# # Used for region-specific counts (not used in current version)
# echo -e "#chr\tstart\tend\tTBR_ID\tTBR_count\tL_count\tR_count\tTBR_specific\tL_specific\tR_specific" > ${OUTFILE}
echo -e "#chr\tstart\tend\tTBR_ID\tTBR_count\tL_count\tR_count" > ${OUTFILE}
paste <( fgrep -v "#" ${TBRs} | cut -f1-4 ) \
<( awk '{ print $NF }' ${COUNTS} | paste - - - | awk -v OFS="\t" '{ print $2, $1, $3 }' ) >> ${OUTFILE}

#Gzip OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
if [ GZI == 1 ]; then
  rm -rf ${TMPI} ${TBRs}
fi
rm ${INTS} ${COUNTS}
