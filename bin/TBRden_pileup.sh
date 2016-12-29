#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Runs intersection of a CNV dataset versus a bed file of elements/loci

#Usage statement
usage(){
cat <<EOF
usage: TBRden_pileup.sh [-h] [-z] [-d distance] [-o OUTFILE] CNVs elements

Runs intersection of a CNV dataset versus a bed file of elements/loci

Positional arguments:
  CNVs      path to CNV input file. Must have at least three columns: chr, 
            start, end
  elements  path to elements input file. Must have at least four columns: chr, 
            element min, element max, and element ID.

Optional arguments:
  -h  HELP          Show this help message and exit
  -d  DIST          Distance padded between window and left/right flanking
                    windows, in bp (default: 1,000,000 bp)
  -z  GZIP          Gzip output file
  -o  OUTFILE       Output file (default: stdout)
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
GZ=0
DIST=1000000
while getopts ":o:m:d:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    d)
      DIST=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
CNVs=$1
elements=$2

#Check for required input
if [ -z ${CNVs} ] || [ -z ${elements} ]; then
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

#Unzip input elements file if gzipped
GZI=0
if [ $( file ${elements} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI=1
  elementsI=`mktemp`; mv ${elementsI} ${elementsI}.gz; elementsI=${elementsI}.gz
  cp ${elements} ${elementsI}
  gunzip ${elementsI}
  elements=$( echo "${elementsI}" | sed 's/\.gz/\t/g' | cut -f1 )
fi

#Make master list of all intervals needed to do pileup calculations
INTS=`mktemp`
fgrep -v "#" ${elements} | awk -v OFS="\t" -v d=${DIST} \
'{ print $1, $2-d, $3-d, $4"_input_"NR"_A\n"$1, $2, $3, $4"_input_"NR"_B\n"$1, $2+d, $3+d, $4"_input_"NR"_C" }' | \
awk -v OFS="\t" -v d=${DIST} '{ if ($2<0) $2=0; if ($3<d) $3=d; print }' > ${INTS}

#Run overlap of master interval file vs CNVs
COUNTS=`mktemp`
bedtools intersect -c -a ${INTS} -b ${CNVs} > ${COUNTS}

#Write results to OUTFILE with header
if [ $( awk '{ print NF }' ${elements} | head -n1 ) -gt 3 ]; then
  echo -e "#chr\tstart\tend\telement_ID\telement_count\tL_count\tR_count" > ${OUTFILE}
  paste <( fgrep -v "#" ${elements} | cut -f1-4 ) \
  <( awk '{ print $NF }' ${COUNTS} | paste - - - | awk -v OFS="\t" '{ print $2, $1, $3 }' ) >> ${OUTFILE}
else
  echo -e "#chr\tstart\tend\telement_count\tL_count\tR_count" > ${OUTFILE}
  paste <( fgrep -v "#" ${elements} | cut -f1-4 ) \
  <( awk '{ print $NF }' ${COUNTS} | paste - - - | awk -v OFS="\t" '{ print $2, $1, $3 }' ) >> ${OUTFILE}
fi

#Gzip OUTFILE, if optioned
if [ ${GZ}==1 ] && [ ${OUTFILE} != "/dev/stdout" ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
if [ GZI == 1 ]; then
  rm -rf ${TMPI} ${TBRs}
fi
rm ${INTS} ${COUNTS}
