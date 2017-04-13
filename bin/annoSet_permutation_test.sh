#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Performs a permutation-based burden test for case vs control CNVs at a \
# set of genomic annotations

#Usage statement
usage(){
cat <<EOF
usage: annoSet_permutation_test.sh [-h] [-N TIMES] [-x EXCLUDE] [-o OUTFILE] CONTROLS CASES ANNO GENOME

Permutation test of CNV burden at a set of genomic annotations by annotation shuffling

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  ANNO       path to input file of annotations to test. Should be BED-style; 
             i.e. must have at least three columns: chr, start, end
  GENOME     BEDTools-style genome file (tab delimmed two-column file: chr, length)

Optional arguments:
  -h  HELP          Show this help message and exit
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -x  EXCLUDE       BED-style intervals to exclude when simulating intervals
  -o  OUTFILE       Output file (default: /dev/stdout)
EOF
}

#Parse arguments
TIMES=1000
OUTFILE=/dev/stdout
EXCLUDE=0
while getopts ":N:x:o:h" opt; do
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
    o)
      OUTFILE=${OPTARG}
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
ANNO=$3
GENOME=$4

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${ANNO} ] || [ -z ${GENOME} ]; then
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

#Unzip input annotation set file if gzipped
GZI_ANNO=0
if [ $( file ${ANNO} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_ANNO=1
  ANNOS=`mktemp`; mv ${ANNOS} ${ANNOS}.gz; ANNOS=${ANNOS}.gz
  cp ${ANNO} ${ANNOS}
  gunzip ${ANNOS}
  ANNOS=$( echo "${ANNOS}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  ANNOS=${ANNO}
fi

#Restrict annotation set to autosomes
grep -e '^[0-9]\|^chr[0-9]' ${ANNOS} > ${ANNOS}2; mv ${ANNOS}2 ${ANNOS}

#Get baseline dCNV
baseline=$( paste <( bedtools intersect -c -a ${ANNOS} -b ${CASE} | \
    awk -v OFS="\t" '{ print $1, $2, $3, $NF }' ) \
  <( bedtools intersect -c -a ${ANNOS} -b ${CTRL} | awk '{ print $NF }' ) | \
  awk -v OFS="\t" '{ sum+=$4-$5 } END { print sum }' )

#Make temporary files for iterating permutations
ANNOS_SHUF=`mktemp`
PERM_OUTPUT=`mktemp`

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  #Print status
  echo "Beginning permutation ${i} of ${TIMES}"

  #Shuffle annotations
  if [ ${EXCLUDE} != "0" ]; then
    bedtools shuffle -f 0.1 -excl ${EXCLUDE} \
    -i ${ANNOS} -g ${GENOME} > ${ANNOS_SHUF}
  else
    bedtools shuffle -noOverlapping -maxTries 10000 -i ${ANNOS} -g ${GENOME} > ${ANNOS_SHUF}
  fi

  #Count pileup of case/control at each element
  paste <( bedtools intersect -c -a ${ANNOS_SHUF} -b ${CASE} | \
    awk -v OFS="\t" '{ print $1, $2, $3, $NF }' ) \
  <( bedtools intersect -c -a ${ANNOS_SHUF} -b ${CTRL} | awk '{ print $NF }' ) | \
  awk -v OFS="\t" '{ sum+=$4-$5 } END { print sum }' >> ${PERM_OUTPUT}
done

#Parameterize null distribution and compute p-value
RES_STAT=`mktemp`
Rscript -e  "dat <- read.table(\"${PERM_OUTPUT}\",header=F)[,1];\
             mu <- mean(dat); sd <- sd(dat); z <- (${baseline}-mu)/sd; p <- 1-pnorm(z);\
             write.table(data.frame(round(mu,4),round(sd,4),round(z,4),p),\
              \"${RES_STAT}\",col.names=F,row.names=F,quote=F,sep=\"\t\")"

#Print results to outfile
for dummy in 1; do
  echo -e "#observed\tperms_greater\tperms_less_or_equal\texpected_mean\texpected_sd\tdifference\tfold_enrichment\tfold_min_95CI\tfold_max_95CI\tZscore\tpvalue"
  for second in 2; do
    echo ${baseline}
    awk -v baseline=${baseline} '{ if ($1>baseline) print $0 }' ${PERM_OUTPUT} | wc -l
    awk -v baseline=${baseline} '{ if ($1<=baseline) print $0 }' ${PERM_OUTPUT} | wc -l
    cut -f1-2 ${RES_STAT}
    awk -v baseline=${baseline} '{ print baseline-$1 }' ${RES_STAT}
    awk -v baseline=${baseline} '{ print baseline/$1 }' ${RES_STAT}
    awk -v baseline=${baseline} '{ print (baseline/$1)-(1.96*($2/$1)) }' ${RES_STAT}
    awk -v baseline=${baseline} '{ print (baseline/$1)+(1.96*($2/$1)) }' ${RES_STAT}
    cut -f3-4 ${RES_STAT}
  done | paste -s
done > ${OUTFILE}

#Clean up
rm -rf ${TMPDIR}
