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
usage: annoSet_permutation_test.sh [-h] [-q QUIET] [-N TIMES] [-x EXCLUDE] [-L LABEL]
                                   [-o OUTFILE] CONTROLS CASES ANNO GENOME

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
  -q  QUIET         Suppresses (some) standard output
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -x  EXCLUDE       BED-style intervals to exclude when simulating intervals
  -L  LABEL         Label of annotation set, printed to output file (default: "Results")
  -o  OUTFILE       Output file (default: /dev/stdout)
EOF
}

#Parse arguments
QUIET=0
TIMES=1000
EXCLUDE=0
LABEL="Results"
OUTFILE=/dev/stdout
while getopts ":N:x:L:o:hq" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    q)
      QUIET=1
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    x)
      EXCLUDE=${OPTARG}
      ;;
    L)
      LABEL=${OPTARG}
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
CTRL=`mktemp`
GZI_CTRL=0
if [ $( file ${CONTROLS} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_CTRL=1
  mv ${CTRL} ${CTRL}.gz; CTRL=${CTRL}.gz
  cp ${CONTROLS} ${CTRL}
  gunzip ${CTRL}
  CTRL=$( echo "${CTRL}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  cp ${CONTROLS} ${CTRL}
fi

#Unzip input case CNV file if gzipped
CASE=`mktemp`
GZI_CASE=0
if [ $( file ${CASES} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_CASE=1
  mv ${CASE} ${CASE}.gz; CASE=${CASE}.gz
  cp ${CASES} ${CASE}
  gunzip ${CASE}
  CASE=$( echo "${CASE}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  cp ${CASES} ${CASE}
fi

#Unzip input annotation set file if gzipped
ANNOS=`mktemp`
GZI_ANNO=0
if [ $( file ${ANNO} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI_ANNO=1
  mv ${ANNOS} ${ANNOS}.gz; ANNOS=${ANNOS}.gz
  cp ${ANNO} ${ANNOS}
  gunzip ${ANNOS}
  ANNOS=$( echo "${ANNOS}" | sed 's/\.gz/\t/g' | cut -f1 )
else
  cp ${ANNO} ${ANNOS}
fi

#Restrict annotation set to autosomes
grep -e '^[0-9]\|^chr[0-9]' ${ANNOS} > ${ANNOS}2
mv ${ANNOS}2 ${ANNOS}

#Get baseline dCNV and case/control counts
base_case=$( bedtools intersect -c -a ${ANNOS} -b ${CASE} | \
  awk -v OFS="\t" '{ print $1, $2, $3, $NF }' | awk '{ sum+=$4 }END{ print sum }' )
base_ctrl=$( bedtools intersect -c -a ${ANNOS} -b ${CTRL} | \
  awk -v OFS="\t" '{ print $1, $2, $3, $NF }' | awk '{ sum+=$4 }END{ print sum }' )
baseline=$( echo -e "${base_case}\t${base_ctrl}" | awk -v OFS="\t" '{ print $1-$2 }' )

#Make temporary files for iterating permutations
ANNOS_SHUF=`mktemp`
PERM_OUTPUT=`mktemp`

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  if [ ${QUIET} == 0 ]; then
    #Print status
    echo "Beginning permutation ${i} of ${TIMES}"
  fi

  #Shuffle annotations
  if [ ${EXCLUDE} != "0" ]; then
    bedtools shuffle -noOverlapping -maxTries 1000 -f 0.1 -excl ${EXCLUDE} \
    -i ${ANNOS} -g ${GENOME} > ${ANNOS_SHUF}
    awk -v OFS="\t" '{ if ($3<$2) $2=$3; print }' ${ANNOS_SHUF} > ${ANNOS_SHUF}2
    mv ${ANNOS_SHUF}2 ${ANNOS_SHUF}
  else
    bedtools shuffle -noOverlapping -maxTries 1000 -i ${ANNOS} -g ${GENOME} > ${ANNOS_SHUF}
    awk -v OFS="\t" '{ if ($3<$2) $2=$3; print }' ${ANNOS_SHUF} > ${ANNOS_SHUF}2
    mv ${ANNOS_SHUF}2 ${ANNOS_SHUF}
  fi

  #Count pileup of case/control at each element
  perm_case=$( bedtools intersect -c -a ${ANNOS_SHUF} -b ${CASE} | \
    awk -v OFS="\t" '{ print $1, $2, $3, $NF }' | awk '{ sum+=$4 }END{ print sum }' )
  perm_ctrl=$( bedtools intersect -c -a ${ANNOS_SHUF} -b ${CTRL} | \
    awk -v OFS="\t" '{ print $1, $2, $3, $NF }' | awk '{ sum+=$4 }END{ print sum }' )
  perm_dCNV=$( echo -e "${perm_case}\t${perm_ctrl}" | awk -v OFS="\t" '{ print $1-$2 }' )
  echo -e "${perm_case}\t${perm_ctrl}\t${perm_dCNV}" >> ${PERM_OUTPUT}
done

#Parameterize null distribution and compute permutation stats
RES_STAT=`mktemp`
#Output columns:
# 1: observed case CNVs
# 2: observed control CNVs
# 3: expected case CNVs
# 4: expected control CNVs
# 5: standard deviation of expected case CNVs
# 6: standard deviation of expected control CNVs
# 7: observed case:control ratio
# 8: expected case:control ratio
# 9: standard deviation of expected case:control ratio
# 10: fold-change of observed vs expected
# 11: lower 95% confidence interval of obs vs exp fold-change
# 12: upper 95% confidence interval of obs vs exp fold-change
# 13: observed Z-score
# 14: p-value (uncorrected)
Rscript -e  "dat <- read.table(\"${PERM_OUTPUT}\",header=F); dat[,4] <- dat[,1]/dat[,2];\
             obs.fold <- ${base_case}/${base_ctrl}; mu <- apply(dat,2,mean); sd <- apply(dat,2,sd);\
             z <- (obs.fold-mu[4])/sd[4]; p <- 1-pnorm(z); fold.est <- obs.fold/mu[4]; \
             fold.lower <- fold.est-(1.96*sd[4]); fold.upper <- fold.est+(1.96*sd[4]); \
             write.table(data.frame(${base_case},${base_ctrl},mu[1],mu[2],sd[1],sd[2],\
             obs.fold,mu[4],sd[4],fold.est,fold.lower,fold.upper,z,p),\
              \"${RES_STAT}\",col.names=F,row.names=F,quote=F,sep=\"\t\")"

#Print results to outfile
for dummy in 1; do
  echo -e "#test\tcase_observed\tcontrol_observed\tcase_expected\tcontrol_expected\t\
case_expected_sd\tcontrol_expected_sd\tobs_case_vs_control\texp_case_vs_control\t\
exp_case_cs_control_sd\tobs_vs_exp\tlower_CI\tupper_CI\tZscore\tp"
  paste <( echo "${LABEL}" ) ${RES_STAT}
done > ${OUTFILE}

#Clean up
rm -rf ${TMPDIR}
