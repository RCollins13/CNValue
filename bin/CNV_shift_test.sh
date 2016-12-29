#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Permutation test of CNV burdens by local CNV shifting

#Usage statement
usage(){
cat <<EOF
usage: CNV_shift_test.sh [-h] [-d DIST] [-N times] [-o OUTDIR]
                         [-p prefix] [-z] CONTROLS CASES BINS

Permutation test of CNV burdens by local CNV shifting

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  BINS       path to input file of bins to test. Should resemble output of TBRden
             burden test output; i.e. must have at least ten columns: chr, start
             coordinate, end coordinate, unique bin ID, CNVs in bin in controls,
             CNVs in left flank in controls, CNVs in right flank in controls, 
             CNVs in bin in cases, CNVs in left flank in controls, and CNVs in
             right flank in controls

Optional arguments:
  -h  HELP          Show this help message and exit
  -d  DIST          Distance to shift each CNV, as fraction of CNV size (default: 2.0)
  -b  BUFFER        Distance padded between window and left/right flanking
                    windows, in bp (default: 1,000,000 bp)
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -o  OUTFILE       Output file (default: stdout)
  -p  PREFIX        Prefix for all output files (default: TBRden_shuffle)
  -z  GZIP          Gzip output
EOF
}

#Parse arguments
DIST=2
BUFFER=1000000
TIMES=1000
OUTFILE=/dev/stdout
PREFIX="TBRden_CNV_shiftTest"
GZ=0
while getopts ":d:b:N:o:p:zh" opt; do
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
    o)
      OUTFILE=${OPTARG}
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
  Rscript -e "write.table(as.data.frame(sample(c(-${DIST},${DIST}),${nCASE},replace=T)),\
              \"${DIRECTION}\",row.names=F,col.names=F,quote=F)"
  paste <( fgrep -v "#" ${CASE} ) ${DIRECTION} | awk -v OFS="\t" \
  '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
  awk -v OFS="\t" '{ if ($2>=0) print }' > ${CASE_SHUF}

  #Shift control CNVs by half their size
  Rscript -e "write.table(as.data.frame(sample(c(-${DIST},${DIST}),${nCTRL},replace=T)),\
              \"${DIRECTION}\",row.names=F,col.names=F,quote=F)"
  paste <( fgrep -v "#" ${CTRL} ) ${DIRECTION} | awk -v OFS="\t" \
  '{ printf "%s\t%i\t%i\n", $1, $2+($NF*($3-$2)), $3+($NF*($3-$2)) }' | \
  awk -v OFS="\t" '{ if ($2>=0) print }' > ${CTRL_SHUF}

  #Run CNV pileups on shuffled CNVs
  ${TBRden_bin}/TBRden_pileup.sh -d ${BUFFER} -o ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CASE_CNVs.pileup.bed ${CASE_SHUF} ${BIN}
  ${TBRden_bin}/TBRden_pileup.sh -d ${BUFFER} -o ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CONTROL_CNVs.pileup.bed ${CTRL_SHUF} ${BIN}

  #Calculate burden from shuffled CNVs
  ${TBRden_bin}/TBRden_test_noPlots.R \
  ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CONTROL_CNVs.pileup.bed.gz \
  ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}.CASE_CNVs.pileup.bed.gz \
  ${TMPDIR}/${PREFIX}_${i} ${PREFIX}_${i}
done

#Collect results
P_RES=`mktemp`
for i in $( seq 1 ${TIMES} ); do
  sed '1d' ${TMPDIR}/${PREFIX}_${i}/${PREFIX}_${i}.TBRden_results.bed | \
  awk '{ print $NF }' | paste -s
done > ${P_RES}

#Compare results with original p-values & print to outfile
echo -e "#chr\tstart\tend\tbin_ID\toriginal_p\tperms_less_sig\tperms_as_or_more_sig\tperm_p" > ${OUTFILE}
for i in $( seq 1 $( cat ${BINS} | wc -l ) ); do
  original=$( sed -n "${i}p" ${BIN} | awk '{ print $NF }' )
  less_sig=$( awk -v i=${i} -v original=${original} \
    '{ if ($(i)>original) print $0 }' ${P_RES} | wc -l )
  as_more_sig=$( awk -v i=${i} -v original=${original} \
    '{ if ($(i)<=original) print $0 }' ${P_RES} | wc -l )
  sed -n "${i}p" ${BINS} | awk -v OFS="\t" -v less=${less_sig} -v as_more=${as_more_sig} \
  '{ print $1, $2, $3, $4, $NF, less, as_more, (1+as_more)/(1+less+as_more) }' 
done >> ${OUTFILE}

#Gzip, if optioned
if [ ${GZ} -eq 1 ]; then
  gzip ${OUTFILE}
fi

#Clean up
rm -rf ${TMPDIR}
