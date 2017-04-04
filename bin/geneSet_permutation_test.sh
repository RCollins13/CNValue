#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Performs a permutation-based burden test for case vs control CNVs for a gene set

#Usage statement
usage(){
cat <<EOF
usage: geneSet_permutation_test.sh [-h] [-N TIMES] [-U UNIVERSE] [-W WHOLE GENE]
                                   [-A ALLOSOMES] [-o OUTFILE] CONTROLS CASES GENESET GTF

Permutation test of CNV burden at a set of genes, normalized by gene length & exonic bases

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  GENESET    path to gene set input file to be tested. Must have one column 
             with gene symbols
  GTF        path to GTF input file. Exons or gene boundaries will be extracted
             from this file

Optional arguments:
  -h  HELP          Show this help message and exit
  -N  TIMES         Number of permutations to perform (default: 1,000)
  -U  UNIVERSE      List of gene symbols to consider as the background 
                    geneset (default: all genes)
  -W  WHOLE GENE    Restrict analysis to CNVs that span the entire gene
                    (default: count any exonic overlap)
  -A  ALLOSOMES     Include allosomes in analyses (default: false)
  -H  OVERRIDE      Path to correctly formatted exons & boundaries directory
                    Note: do not use this option unless you know what you're doing
  -o  OUTFILE       Output file (default: /dev/stdout)
EOF
}

#Parse arguments
TIMES=1000
OUTFILE=/dev/stdout
UNIVERSE=ALL
WG=0
ALLO=0
OVER=0
while getopts ":N:U:WA:H:o:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    N)
      TIMES=${OPTARG}
      ;;
    U)
      UNIVERSE=${OPTARG}
      ;;
    W)
      WG=1
      ;;
    A)
      ALLO=1
      ;;
    H)
      OVER=${OPTARG}
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
GENESET=$3
GTF=$4

####NOTE: REPLACES ALL HYPHENS WITH UNDERSCORES IN GENE SYMBOLS FOR GREP COMPATIBILITY####

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${GENESET} ] || [ -z ${GTF} ]; then
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

#Restricts CNVs to autosomes unless optioned
if [ ${ALLO} -eq 0 ]; then
  grep -e '^[0-9]\|^chr[0-9]' ${CTRL} > ${CTRL}2; mv ${CTRL}2 ${CTRL}
  grep -e '^[0-9]\|^chr[0-9]' ${CASE} > ${CASE}2; mv ${CASE}2 ${CASE}
fi

#Skip dictionary creation if hard override is optioned (saves time)
if [ ${OVER} == "0" ]; then
  #Parse exon definitions from GTF
  echo -e "STATUS::$(date)::BUILDING GENE UNIVERSE DICTIONARY FROM GTF..."
  EXONS=`mktemp`
  fgrep -v "#" ${GTF} | sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
  '{ if ($3=="exon") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
  awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
  sed 's/^chr//g' | sed 's/\-/_/g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > ${EXONS}

  #Parse gene boundary definitions from GTF
  BOUNDARIES=`mktemp`
  fgrep -v "#" ${GTF} | sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
  '{ if ($3=="gene") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
  awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
  sed 's/^chr//g' | sed 's/\-/_/g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > ${BOUNDARIES}
else
  EXONS=`mktemp`
  cp ${OVER}/exons.bed ${EXONS}
  BOUNDARIES=`mktemp`
  cp ${OVER}/boundaries.bed ${BOUNDARIES}
fi

#Subset exons and boundaries to autosomes unless optioned
if [ ${ALLO} -eq 0 ]; then
  grep -e '^[0-9]\|^chr[0-9]' ${EXONS} > ${EXONS}2; mv ${EXONS}2 ${EXONS}
  grep -e '^[0-9]\|^chr[0-9]' ${BOUNDARIES} > ${BOUNDARIES}2; mv ${BOUNDARIES}2 ${BOUNDARIES}
fi

#Restrict universal set to gene symbols that appear in exons and/or boundaries
#If no universal set provided, define universal set as list of unique gene 
#symbols with at least one exon and/or gene boundary defined
UNIV=`mktemp`
if [ ${UNIVERSE} != "ALL" ]; then
  sed 's/\-/_/g' ${UNIVERSE} | fgrep -wf - <( cat ${EXONS} ${BOUNDARIES} ) | \
  cut -f4 | sort | uniq > ${UNIV}
else
  cat ${EXONS} ${BOUNDARIES} | cut -f4 | sort | uniq > ${UNIV}
fi

#Subset exon & boundary definitions to genes that exist in the universal set
fgrep -wf ${UNIV} ${EXONS} > ${EXONS}2; mv ${EXONS}2 ${EXONS}
fgrep -wf ${UNIV} ${BOUNDARIES} > ${BOUNDARIES}2; mv ${BOUNDARIES}2 ${BOUNDARIES}

#Substitute dashes for underscores in supplied gene list
QUERY=`mktemp`
sed 's/\-/_/g' ${GENESET} > ${QUERY}

#Count number of query genes with a match in either exon/boundary file
QUERY_IN_UNIV=$( fgrep -wf ${QUERY} ${UNIV} | wc -l )

# #Compute boundary size, exonic bases, and baseline dCNVs per gene
# echo -e "STATUS::$(date)::BUILDING GENE UNIVERSE REFERENCE"
# BASELINE=`mktemp`
# while read gene; do
#   for dummy in 1; do
#     echo ${gene}
#     #Boundary size
#     awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${BOUNDARIES} | \
#     sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
#     awk '{ sum+=$3-$2 }END{ print sum }'
#     #Exonic bases
#     awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${EXONS} | \
#     sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
#     awk '{ sum+=$3-$2 }END{ print sum }'
#     if [ ${WG} -eq 1 ]; then
#       #Boundary-based dCNV
#       paste <( awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${BOUNDARIES} | \
#                bedtools intersect -wb -f 1 -a - -b ${CASE} | cut -f5- | \
#                sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) \
#             <( awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${BOUNDARIES} | \
#                bedtools intersect -wb -f 1 -a - -b ${CTRL} | cut -f5- | \
#                sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) | \
#       awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }'
#     else
#       #Exon-based dCNV
#       paste <( awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${EXONS} | \
#                bedtools intersect -u -a ${CASE} -b - | wc -l ) \
#             <( awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${EXONS} | \
#                bedtools intersect -u -a ${CTRL} -b - | wc -l ) | \
#       awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }'
#     fi
#   done | paste -s 
# done < ${UNIVERSE} > ${BASELINE}

#Compute baseline dCNV for query set
if [ ${WG} -eq 1 ]; then
  #Boundary-based dCNV
  OBS_dCNV=$( paste <( fgrep -wf ${QUERY} ${BOUNDARIES} | \
           bedtools intersect -wb -f 1 -a - -b ${CASE} | cut -f5- | \
           sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) \
        <( fgrep -wf ${QUERY} ${BOUNDARIES} | \
           bedtools intersect -wb -f 1 -a - -b ${CTRL} | cut -f5- | \
           sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) | \
  awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }' )
else
  #Exon-based dCNV
  OBS_dCNV=$( paste <( fgrep -wf ${QUERY} ${EXONS} | \
           bedtools intersect -u -a ${CASE} -b - | wc -l ) \
        <( fgrep -wf ${QUERY} ${EXONS} | \
           bedtools intersect -u -a ${CTRL} -b - | wc -l ) | \
  awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }' )
fi

#Make temporary files for iterating permutations
SHUF_GENES=`mktemp`
PERM_OUTPUT=`mktemp`

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  #Print status
  echo "Beginning permutation ${i} of ${TIMES}"

  #Sample new set of genes
  shuf ${UNIV} | head -n${QUERY_IN_UNIV} > ${SHUF_GENES}

  #Recompute dCNV
  if [ ${WG} -eq 1 ]; then
    #Boundary-based dCNV
    paste <( fgrep -wf ${SHUF_GENES} ${BOUNDARIES} | \
             bedtools intersect -wb -f 1 -a - -b ${CASE} | cut -f5- | \
             sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) \
          <( fgrep -wf ${SHUF_GENES} ${BOUNDARIES} | \
             bedtools intersect -wb -f 1 -a - -b ${CTRL} | cut -f5- | \
             sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq | wc -l ) | \
    awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }'
  else
    #Exon-based dCNV
    paste <( fgrep -wf ${SHUF_GENES} ${EXONS} | \
             bedtools intersect -u -a ${CASE} -b - | wc -l ) \
          <( fgrep -wf ${SHUF_GENES} ${EXONS} | \
             bedtools intersect -u -a ${CTRL} -b - | wc -l ) | \
    awk -v OFS="\t" '{ sum+=$1-$2 } END { print sum }'
  fi >> ${PERM_OUTPUT}
done

#Parameterize null distribution and compute p-value
RES_STAT=`mktemp`
Rscript -e  "dat <- read.table(\"${PERM_OUTPUT}\",header=F)[,1];\
             mu <- mean(dat); sd <- sd(dat); z <- (${OBS_dCNV}-mu)/sd; p <- 1-pnorm(z);\
             write.table(data.frame(round(mu,4),round(sd,4),round(z,4),p),\
              \"${RES_STAT}\",col.names=F,row.names=F,quote=F,sep=\"\t\")"

#Print results to outfile
for dummy in 1; do
  echo -e "#observed\tperms_greater\tperms_less_or_equal\texpected_mean\texpected_sd\tdifference\tfold_enrichment\tfold_min_95CI\tfold_max_95CI\tZscore\tpvalue"
  for second in 2; do
    echo ${OBS_dCNV}
    awk -v baseline=${OBS_dCNV} '{ if ($1>baseline) print $0 }' ${PERM_OUTPUT} | wc -l
    awk -v baseline=${OBS_dCNV} '{ if ($1<=baseline) print $0 }' ${PERM_OUTPUT} | wc -l
    cut -f1-2 ${RES_STAT}
    awk -v baseline=${OBS_dCNV} '{ print baseline-$1 }' ${RES_STAT}
    awk -v baseline=${OBS_dCNV} '{ print baseline/$1 }' ${RES_STAT}
    awk -v baseline=${OBS_dCNV} '{ print (baseline/$1)-(1.96*($2/$1)) }' ${RES_STAT}
    awk -v baseline=${OBS_dCNV} '{ print (baseline/$1)+(1.96*($2/$1)) }' ${RES_STAT}
    cut -f3-4 ${RES_STAT}
  done | paste -s
done > ${OUTFILE}

#Clean up
rm -rf ${TMPDIR}
