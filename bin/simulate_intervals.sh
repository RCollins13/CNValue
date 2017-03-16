#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Permutation test of CNV burdens by local CNV shifting

#Usage statement
usage(){
cat <<EOF
usage: simulate_intervals.sh [-h] [-s SIZE] [-d STDEV] [-N COUNT] [-x EXCLUDE] [-o OUTFILE] [-z] GENOME

Simulates BED-style intervals whose sizes are normally distributed

Positional arguments:
  GENOME     BEDTools-style genome file (tab delimmed two-column file: chr, length)

Optional arguments:
  -h  HELP          Show this help message and exit
  -s  SIZE          Average size of simulated intervals (default: 5000bp)
  -d  STDEV         Standard deviation of simulated interval sizes (default: 1000bp)
  -N  COUNT         Number of intervals to simulate (default: 10000)
  -x  EXCLUDE       BED-style intervals to exclude when simulating intervals
  -o  OUTFILE       Output file (default: stdout)
  -z  GZIP          Gzip output
EOF
}

#Parse arguments
SIZE=5000
STDEV=1000
COUNT=10000
EXCLUDE=0
OUTFILE=/dev/stdout
GZ=0
while getopts ":s:d:N:x:o:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    s)
      SIZE=${OPTARG}
      ;;
    d)
      STDEV=${OPTARG}
      ;;
    N)
      COUNT=${OPTARG}
      ;;
    x)
      EXCLUDE=${OPTARG}
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
GENOME=$1

#Check that genome file exists
if [ -s ${GENOME} ]; then
  echo -e "\nERROR: GENOME FILE DOES NOT EXIST\n"
  usage
  exit 0
fi

#Subtract exclusion intervals from genome file if specified
SAMPSPACE=`mktemp`
if [ ${EXCLUDE} != 0 ]; then
  if ! [ -s ${EXCLUDE} ]; then
    echo -e "\nERROR: -x IS SPECIFIED BUT EXCLUSION FILE DOES NOT EXIST\n"
    usage
    exit 0
  else
    bedtools subtract -a <( awk -v OFS="\t" '{ print $1, "1", $2 }' ${GENOME} ) \
    -b ${EXCLUDE} > ${SAMPSPACE}
  fi
else
  awk -v OFS="\t" '{ print $1, "1", $2 }' ${GENOME} > ${SAMPSPACE}
fi
NGENOME=`mktemp`

#Compute length modifiers for intervals (5x total count requested)
MODS=`mktemp`
Rscript -e "write.table(round(rnorm(n=5*${COUNT},sd=${STDEV})),\"${MODS}\",\
  col.names=F,row.names=F,sep=\"\t\",quote=F)"

#Generate random intervals
bedtools random -l ${SIZE} -n $((5*${COUNT})) -g ${GENOME} | paste - ${MODS} | \
awk -v OFS="\t" '{ print $1, $2+$NF, $3 }' | awk -v OFS="\t" \
'{ if ($3-$2>=1 && $2>0) print $0 }' | bedtools intersect -f 1 -a - -b ${SAMPSPACE} | \
head -n ${COUNT} | sort -Vk1,1 -k2,2n -k3,3n > ${OUTFILE}

#Gzip output if optioned
if [ ${OUTFILE} != "/dev/stdout" ] && [ ${GZ} -gt 0 ]; then
  gzip -f ${OUTFILE}
fi

#Clean up
rm -rf ${SAMPSPACE} ${MODS}
