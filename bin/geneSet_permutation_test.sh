#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Performs a permutation-based burden test for case vs control CNVs for a gene set

#Testing dev parameters
TIMES=100
OUTFILE=${TMPDIR}/geneset_test.out
UNIVERSE=${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list
WG=0
ALLO=0
OVER=${WRKDIR}/data/misc/exons_boundaries_dictionary/
CONTROLS=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.E3.GRCh37.all.bed.gz
CASES=${WRKDIR}/data/CNV/CNV_MASTER/NDD/NDD.DEL.E3.GRCh37.all.bed.gz
GENESET=${WRKDIR}/data/master_annotations/genelists/DDD_2017.genes.list
GTF=${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf
QUIET=0
LABEL="DDD_2017"

#Usage statement
usage(){
cat <<EOF
usage: geneSet_permutation_test.sh [-h] [-N TIMES] [-U UNIVERSE] [-W WHOLE GENE]
                                   [-A ALLOSOMES] [-H OVERRIDE] [-o OUTFILE] [-q QUIET] 
                                   CONTROLS CASES GENESET GTF

Permutation test of CNV burden at a set of genes, normalized by gene length & exonic bases

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  GENESET    path to gene set input file to be tested. Must have one column 
             with gene symbols
  GTF        path to GTF input file. Exons or gene boundaries will be extracted
             from this file. Can be overridden with -H

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
  -L  LABEL         Label of annotation set, printed to output file (default: "Results")
  -o  OUTFILE       Output file (default: /dev/stdout)
  -q  QUIET         Suppresses (some) standard output
EOF
}

#Parse arguments
TIMES=1000
LABEL="Results"
OUTFILE=/dev/stdout
UNIVERSE=ALL
WG=0
ALLO=0
OVER=0
QUIET=0
while getopts ":N:U:WA:H:L:o:hq" opt; do
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
    L)
      LABEL=${OPTARG}
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
    q)
      QUIET=1
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
  if [ ${QUIET} -eq 0 ]; then
    echo -e "STATUS::$(date)::BUILDING GENE UNIVERSE DICTIONARY FROM GTF..."
  fi
  EXONS=`mktemp`
  fgrep -v "#" ${GTF} | sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
  '{ if ($3=="exon") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
  awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
  sed 's/^chr//g' | sed 's/\-/_/g' | \
  sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > ${EXONS}
  #Parse gene boundary definitions from GTF
  BOUNDARIES=`mktemp`
  fgrep -v "#" ${GTF} | sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
  '{ if ($3=="gene") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
  awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
  sed 's/^chr//g' | sed 's/\-/_/g' | \
  sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > ${BOUNDARIES}
else
  if [ ${QUIET} -eq 0 ]; then
    echo -e "STATUS::$(date)::LOADING GENE UNIVERSE DICTIONARY FROM PATH..."
  fi
  EXONS=`mktemp`
  sed 's/^chr//g' ${OVER}/exons.bed | sed 's/\-/_/g' | \
  sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > ${EXONS}
  BOUNDARIES=`mktemp`
  sed 's/^chr//g' ${OVER}/boundaries.bed | sed 's/\-/_/g' | \
  sort -Vk1,1 -k2,2n -k3,3n -k4,4 | uniq > ${BOUNDARIES}
fi

#Subset exons and boundaries to autosomes unless optioned
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::HARMONIZING INPUT DATASETS..."
fi
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
  awk '{ print $4 }' | sort | uniq > ${UNIV}
else
  cat ${EXONS} ${BOUNDARIES} | awk '{ print $4 }' | sort | uniq > ${UNIV}
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

#Calculate unique case & control CNV counts for all genes in universe
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::COMPUTING BASELINE CNV BURDENS..."
fi
ALL_GENES_CASE_CNV=`mktemp`
ALL_GENES_CTRL_CNV=`mktemp`
if [ ${WG} -eq 1 ]; then
  #Case, whole-gene CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${BOUNDARIES} -b ${CASE} | \
  awk -v OFS="\t" '{ print $4, $8 }' | \
  sort -k1,1 -k2,2 | uniq | cut -f1 | uniq -c | \
  awk -v OFS="\t" '{ print $2, $1 }' > ${ALL_GENES_CASE_CNV}
  #Control, whole-gene CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${BOUNDARIES} -b ${CTRL} | \
  awk -v OFS="\t" '{ print $4, $8 }' | \
  sort -k1,1 -k2,2 | uniq | cut -f1 | uniq -c | \
  awk -v OFS="\t" '{ print $2, $1 }' > ${ALL_GENES_CTRL_CNV}
else
  #Case, exonic CNVs
  bedtools intersect -wa -wb -a ${CASE} -b ${EXONS} | \
  awk -v OFS="\t" '{ print $NF, $4 }' | sort -k1,1 -k2,2 | uniq | \
  cut -f1 | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' > ${ALL_GENES_CASE_CNV}
  #Control, exonic CNVs
  bedtools intersect -wa -wb -a ${CTRL} -b ${EXONS} | \
  awk -v OFS="\t" '{ print $NF, $4 }' | sort -k1,1 -k2,2 | uniq | \
  cut -f1 | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' > ${ALL_GENES_CTRL_CNV}
fi

#Compute baseline dCNV and case/control counts for query set
base_case=$( fgrep -wf ${QUERY} ${ALL_GENES_CASE_CNV} | \
             awk '{ sum+=$2 }END{ print sum }' )
base_ctrl=$( fgrep -wf ${QUERY} ${ALL_GENES_CTRL_CNV} | \
             awk '{ sum+=$2 }END{ print sum }' )
baseline=$( echo -e "${base_case}\t${base_ctrl}" | \
            awk -v OFS="\t" '{ print $1-$2 }' )

#Make temporary files for iterating permutations
GENES_SHUF=`mktemp`
PERM_OUTPUT=`mktemp`

#Repeat for number of permutations specified by user
for i in $( seq 1 ${TIMES} ); do
  if [ ${QUIET} == 0 ]; then
    #Print status
    echo "Beginning permutation ${i} of ${TIMES}"
  fi

  #Sample new set of genes
  shuf ${UNIV} | head -n${QUERY_IN_UNIV} > ${GENES_SHUF}

  #Count pileup of case/control at each element
  perm_case=$( fgrep -wf ${GENES_SHUF} ${ALL_GENES_CASE_CNV} | \
               awk '{ sum+=$2 }END{ print sum }' )
  perm_ctrl=$( fgrep -wf ${GENES_SHUF} ${ALL_GENES_CTRL_CNV} | \
               awk '{ sum+=$2 }END{ print sum }' )
  perm_dCNV=$( echo -e "${perm_case}-${perm_ctrl}" | bc )

  #Print permutation results to file
  echo -e "${perm_case}\t${perm_ctrl}\t${perm_dCNV}" >> ${PERM_OUTPUT}
done

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
exp_case_vs_control_sd\tobs_vs_exp\tlower_CI\tupper_CI\tZscore\tp"
  paste <( echo "${LABEL}" ) ${RES_STAT}
done > ${OUTFILE}

#Clean up
rm -rf ${TMPDIR}
