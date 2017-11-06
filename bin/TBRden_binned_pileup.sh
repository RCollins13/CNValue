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
usage: TBRden_binned_pileup.sh [-h] [-w WINDOW] [-s STEP] [-d DIST] [-r SMOOTH] [-R WINDOW_SMOOTH]
                               [-I OVERLAP] [-x EXCLUDE] [-o OUTFILE] [-z] CNVs genome

Runs intersection of a CNV dataset versus a binned genome of sliding windows

Positional arguments:
  CNVs     path to CNV input file. Must have at least three columns: chr, start, end
  genome   path to input genome file (per bedtools specs). Must have two columns:
           contig ID & length

Optional arguments:
  -h  HELP          Show this help message and exit
  -z  GZIP          Gzip output file
  -w  WINDOW        Size of window, in bp (default: 5,000 bp)
  -s  STEP          Size of step, in bp (default: 5,000 bp)
  -d  DIST          Distance padded between window and left/right flanking
                    windows, in bp (default: 1,000,000 bp)
  -r  SMOOTH        Number of bins to pad to each bin for smoothing (default: 0)
  -R  WINDOW_SMOOTH Number of bins to pad to each window bin for smoothing (default: 0)
  -I  OVERLAP       Minimum overlap of bin required by each CNV, as pct (default: 0)
  -x  EXCLUDE       Regions to exclude from calculations
  -o  OUTFILE       Output file (default: stdout)
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
GZ=0
EXCLUDE=0
WINDOW=5000
STEP=5000
DIST=1000000
SMOOTH=0
WSMOOTH=0
OVERLAP=0
while getopts ":o:w:s:d:r:R:I:x:zh" opt; do
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
    r)
      SMOOTH=${OPTARG}
      ;;
    R)
      WSMOOTH=${OPTARG}
      ;;
    I)
      OVERLAP=${OPTARG}
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

#DEV PARAMETERS
WINDOW=200000
GZIP=1
STEP=10000
DIST=2500000
WSMOOTH=2
SMOOTH=0
OVERLAP=1
pheno=NDD
CNV=DEL
VF=E2
filt=all
OUTFILE=${TMPDIR}/test_out.bed
EXCLUDE=${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed
CNVs=${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz
genome=/data/talkowski/rlc47/src/GRCh37.genome

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
  done < ${genome} > ${INTS}
else
  while read contig length; do
    paste <( seq 0 ${STEP} $(( ${length}-${WINDOW} )) ) \
    <( seq ${WINDOW} ${STEP} ${length} ) | \
    awk -v OFS="\t" -v contig=${contig} '{ print contig, $1, $2 }' | awk -v OFS="\t" \
    '{ print $0, "Bin_"$1"_"NR }'
  done < ${genome} > ${INTS}
fi

#Run overlap of master interval file vs CNVs
COUNTS=`mktemp`
bedtools intersect -c -f ${OVERLAP} -a ${INTS} -b ${CNVs} > ${COUNTS}

#Get left flanking medians
LEFT=`mktemp`
awk -v OFS="\t" -v d=${DIST} -v buffer=$(( ${WSMOOTH}*${WINDOW} )) \
'{ print $1, $2-d-buffer, $3-d+buffer, $4"_L" }' ${INTS} | \
awk -v OFS="\t" -v w=${WINDOW} '{ if ($2<0) $2=0; if ($3<w) $3=w; print }' | \
bedtools map -c 5 -o median -a - -b ${COUNTS} | awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | \
cut -f1 -d\. > ${LEFT}

#Get right flanking medians
RIGHT=`mktemp`
awk -v OFS="\t" -v d=${DIST} -v buffer=$(( ${WSMOOTH}*${WINDOW} )) \
'{ print $1, $2+d-buffer, $3+d+buffer, $4"_R" }' ${INTS} | \
awk -v OFS="\t" -v w=${WINDOW} '{ if ($2<0) $2=0; if ($3<w) $3=w; print }' | \
bedtools map -c 5 -o median -a - -b ${COUNTS} | awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | \
cut -f1 -d\. > ${RIGHT}

#Smooth main bin values if optioned
if [ ${SMOOTH} -gt 0 ]; then
  for file in ${COUNTS} ${LEFT} ${RIGHT}; do
    awk -v S=${SMOOTH} -v w=${WINDOW} -v OFS="\t" '{ print $1, $2-(w*S), $3+(w*S), $4 }' ${file} | \
    awk -v OFS="\t" -v w=${WINDOW} '{ if ($2<0) $2=0; if ($3<w) $3=w; print }' | \
    bedtools map -c 5 -o mean -a - -b <( awk '{ if ($5!="NA") print $0 }' ${file} ) | \
    awk -v S=${SMOOTH} -v w=${WINDOW} -v OFS="\t" '{ print $1, $2+(w*S), $3-(w*S), $4, $5 }' | \
    awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' > ${file}2
    mv ${file}2 ${file}
  done
fi

#Write results to OUTFILE with header
echo -e "#chr\tstart\tend\tBIN_ID\tBIN\tLEFT\tRIGHT" > ${OUTFILE}
paste ${COUNTS} <( cut -f5 ${LEFT} ) <( cut -f5 ${RIGHT} ) >> ${OUTFILE}

#Gzip OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
rm ${INTS} ${COUNTS} ${LEFT} ${RIGHT}
