#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Runs intersection of a CNV dataset versus a binned genome of sliding windows

#Usage statement
usage(){
cat <<EOF
usage: TBRden_binned_pileup.sh [-h] [-w WINDOW] [-s STEP] [-d DIST]
                        [-x EXCLUDE] [-o OUTFILE] [-z] CNVs genome

Runs intersection of a CNV dataset versus a binned genome of sliding windows

Positional arguments:
  CNVs     path to CNV input file. Must have at least three columns: chr, start, end
  genome   path to input genome file (per bedtools specs). Must have two columns:
           contig ID & length

Optional arguments:
  -h  HELP          Show this help message and exit
  -z  GZIP          Gzip output file
  -w  WINDOW        Size of window, in bp (default: 100,000 bp)
  -s  STEP          Size of step, in bp (default: 25,000 bp)
  -d  DIST          Distance padded between window and left/right flanking
                    windows, in bp (default: 1,000,000 bp)
  -x  EXCLUDE       Regions to exclude from calculations
  -o  OUTFILE       Output file (default: stdout)
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
GZ=0
EXCLUDE=0
WINDOW=100000
STEP=25000
DIST=1000000
while getopts ":o:w:s:d:x:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    w)
      WINDOW=${OPTARG}
      ;;
    s)
      STEP=${OPTARG}
      ;;
    d)
      DIST=${OPTARG}
      ;;
    x)
      EXCLUDE=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1 ))
CNVs=$1
genome=$2

#Check for required input
if [ -z ${CNVs} ] || [ -z ${genome} ]; then
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

#1. Make master list of all intervals needed to do pileup calculations
#2. Remove all bins with any overlap versus excluded regions
INTS=`mktemp`
if ! [ ${EXCLUDE} == "0" ]; then
  while read contig length; do
    paste <( seq 0 ${STEP} $(( ${length}-${WINDOW} )) ) \
    <( seq ${WINDOW} ${STEP} ${length} ) | \
    awk -v OFS="\t" -v contig=${contig} '{ print contig, $1, $2 }' | \
    bedtools intersect -v -a - -b ${EXCLUDE} | awk -v OFS="\t" \
    '{ print $0, "Bin_"$1"_"NR }'
  done < ${genome} | awk -v OFS="\t" -v W=${WINDOW} -v d=${DIST} \
  '{ print $1, $2-(W+d), $2-d, $4"_A\n"$1, $2, $3, $4"_B\n"$1, $3+d, $3+(W+d), $4"_C" }' > ${INTS}
else
  while read contig length; do
    paste <( seq 0 ${STEP} $(( ${length}-${WINDOW} )) ) \
    <( seq ${WINDOW} ${STEP} ${length} ) | \
    awk -v OFS="\t" -v contig=${contig} '{ print contig, $1, $2 }' | awk -v OFS="\t" \
    '{ print $0, "Bin_"$1"_"NR }'
  done < ${genome} | awk -v OFS="\t" -v W=${WINDOW} \
  '{ print $1, $2-(W+d), $2-d, $4"_A\n"$1, $2, $3, $4"_B\n"$1, $3+d, $3+(W+d), $4"_C" }' > ${INTS}
fi

#Reset coordinates of bins with negative size or bins extending beyond the ends of chromosomes to 0
awk -v OFS="\t" -v w=${WINDOW} '{ if ($2<0) $2=0; if ($3<w) $3=w; print }' ${INTS} | \
bedtools intersect -u -a - -b <( awk -v OFS="\t" '{ print $1, 1, $2 }' ${genome} ) > ${INTS}2
awk -v OFS="\t" -v w=${WINDOW} '{ if ($2<0) $2=0; if ($3<w) $3=w; print }' ${INTS} | \
bedtools intersect -v -a - -b <( awk -v OFS="\t" '{ print $1, 1, $2 }' ${genome} ) | \
awk -v OFS="\t" -v w=${WINDOW} '{ print $1, 0, w, $4 }' >> ${INTS}2
sort -Vk1,1 -k4,4V ${INTS}2 > ${INTS}

#Run overlap of master interval file vs CNVs
COUNTS=`mktemp`
bedtools intersect -c -a ${INTS} -b ${CNVs} > ${COUNTS}

#Write results to OUTFILE with header
echo -e "#chr\tstart\tend\tBIN_ID\tBIN_count\tL_count\tR_count" > ${OUTFILE}
paste <( fgrep -v "#" ${INTS} | fgrep "_B" | cut -f1-4 | sed 's/_B//g' ) \
<( awk '{ print $NF }' ${COUNTS} | paste - - - | awk -v OFS="\t" '{ print $2, $1, $3 }' ) >> ${OUTFILE}

#Gzip OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
rm ${INTS} ${COUNTS}
