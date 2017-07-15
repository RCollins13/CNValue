#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects data required for modeling per-gene rCNV burden scores

#Testing dev parameters
# export WRKDIR=/data/talkowski/Samples/rCNVmap
# source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh
# OUTFILE=${TMPDIR}/geneset_test.out
# UNIVERSE=${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list
# WG=0
# ALLO=0
# QUIET=0
# OVER=${WRKDIR}/data/misc/exons_boundaries_dictionary/
# CONTROLS=${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.E3.GRCh37.all.bed.gz
# CASES=${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.DEL.E3.GRCh37.all.bed.gz
# GTF=${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf
# REF=${h37}
# BIN=${WRKDIR}/bin/rCNVmap/bin/

#Usage statement
usage(){
cat <<EOF
usage: gather_geneScore_data.sh [-h] [-U UNIVERSE] [-W WHOLE GENE] [-A ALLOSOMES]
                                [-H OVERRIDE] [-o OUTFILE] [-q QUIET] 
                                CONTROLS CASES GTF REF

Collects data required for modeling per-gene rCNV burden scores

Positional arguments:
  CONTROLS   path to control CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  CASES      path to case CNV input file. Must have at least three columns: 
             chr, CNV start, CNV end
  GTF        path to GTF input file. Exons and gene boundaries will be extracted
             from this file. Can be overridden with -H
  REF        path to reference .fa 

Optional arguments:
  -h  HELP          Show this help message and exit
  -U  UNIVERSE      List of gene symbols to consider as the background 
                    geneset (default: all genes)
  -W  WHOLE GENE    Restrict analysis to CNVs that span the entire gene
                    (default: count any exonic overlap)
  -A  ALLOSOMES     Include allosomes in analyses (default: false)
  -H  OVERRIDE      Path to correctly formatted exons & boundaries directory
                    Note: do not use this option unless you know what you're doing
  -o  OUTFILE       Output file (default: /dev/stdout)
  -q  QUIET         Suppresses (some) standard output
EOF
}

#Parse arguments
OUTFILE=/dev/stdout
UNIVERSE=ALL
WG=0
ALLO=0
OVER=0
QUIET=0
while getopts ":U:WAH:o:hq" opt; do
  case "$opt" in
    h)
      usage
      exit 0
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
    q)
      QUIET=1
      ;;
  esac
done
shift $((${OPTIND} - 1))
CONTROLS=$1
CASES=$2
GTF=$3
REF=$4

#Get path to rCNVmap bin
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

####NOTE: REPLACES ALL HYPHENS WITH UNDERSCORES IN GENE SYMBOLS FOR GREP COMPATIBILITY####

#Check for required input
if [ -z ${CONTROLS} ] || [ -z ${CASES} ] || [ -z ${GTF} ] || [ -z ${REF} ]; then
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

#Annotates gene boundaries file with GC content from reference
GC_PROFILES=`mktemp`
bedtools nuc -fi ${REF} -bed ${BOUNDARIES} | fgrep -v "#" | cut -f4,6 | \
sort -k1,1 > ${GC_PROFILES}

#Isolate unique CNV-gene pairs
CNV_GENE_PAIRS_CASE=`mktemp`
CNV_GENE_PAIRS_CTRL=`mktemp`
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::OVERLAPPING CNVS AND GENES..."
fi
if [ ${WG} -eq 1 ]; then
  #Case, whole-gene CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${BOUNDARIES} -b ${CASE} | \
  awk -v OFS="\t" '{ print $8, $4 }' | \
  sort -k1,1 -k2,2 | uniq > ${CNV_GENE_PAIRS_CASE}
  #Control, whole-gene CNVs
  bedtools intersect -f 1.0 -wa -wb -a ${BOUNDARIES} -b ${CTRL} | \
  awk -v OFS="\t" '{ print $8, $4 }' | \
  sort -k1,1 -k2,2 | uniq > ${CNV_GENE_PAIRS_CTRL}
else
  #Case, exonic CNVs
  bedtools intersect -wa -wb -a ${CASE} -b ${EXONS} | \
  awk -v OFS="\t" '{ print $4, $NF }' | sort -k1,1 -k2,2 | \
  uniq > ${CNV_GENE_PAIRS_CASE}
  #Control, exonic CNVs
  bedtools intersect -wa -wb -a ${CTRL} -b ${EXONS} | \
  awk -v OFS="\t" '{ print $4, $NF }' | sort -k1,1 -k2,2 | \
  uniq > ${CNV_GENE_PAIRS_CTRL}
fi

#Calculate baseline dCNV & dGene info
if [ ${QUIET} -eq 0 ]; then
  echo -e "STATUS::$(date)::RUNNING CNV WEIGHTING MODEL..."
fi
${BIN}/gather_geneScore_data.helper.R \
${CNV_GENE_PAIRS_CASE} \
${CNV_GENE_PAIRS_CTRL} \
${TMPDIR}/

#Write header to output file
echo -e "#gene\tgene_length\texonic_bases\tGC\tcase_CNV\tcontrol_CNV\tcase_CNV_weighted\tcontrol_CNV_weighted" > ${OUTFILE}

#Compute boundary size, exonic bases, and GC per gene
while read gene; do
  #Only compute if gene is included in boundaries file
  if [ $( awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${BOUNDARIES} | wc -l ) -gt 0 ]; then 
    for dummy in 1; do
      echo ${gene}
      #Boundary size
      awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${BOUNDARIES} | \
      sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
      awk '{ sum+=$3-$2 }END{ print sum }'
      #Exonic bases
      awk -v gene=${gene} '{ if ($4==gene) print $0 }' ${EXONS} | \
      sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
      awk '{ sum+=$3-$2 }END{ print sum }'
      #GC content (take average in the case of multiple hits--this is rare)
      awk -v gene=${gene} '{ if ($1==gene) print $2 }' ${GC_PROFILES} | \
      awk '{ sum+=$1 }END{ print sum/NR }'
      #Case CNVs
      caseCNV=$( awk -v gene=${gene} '{ if ($1==gene) print $2 }' ${TMPDIR}/CASE.CNVsPerGene.txt )
      if [ -z ${caseCNV} ]; then
        caseCNV=0
      fi
      echo "${caseCNV}"
      #Control CNVs
      ctrlCNV=$( awk -v gene=${gene} '{ if ($1==gene) print $2 }' ${TMPDIR}/CTRL.CNVsPerGene.txt )
      if [ -z ${ctrlCNV} ]; then
        ctrlCNV=0
      fi
      echo "${ctrlCNV}"
      #Case CNVs (weighted)
      caseCNVweighted=$( awk -v gene=${gene} '{ if ($1==gene) print $2 }' ${TMPDIR}/CASE.weightedCNVsPerGene.txt )
      if [ -z ${caseCNVweighted} ]; then
        caseCNVweighted=0
      fi
      echo "${caseCNVweighted}"
      #Control CNVs (weighted)
      ctrlCNVweighted=$( awk -v gene=${gene} '{ if ($1==gene) print $2 }' ${TMPDIR}/CTRL.weightedCNVsPerGene.txt )
      if [ -z ${ctrlCNVweighted} ]; then
        ctrlCNVweighted=0
      fi
      echo "${ctrlCNVweighted}"
    done | paste -s
  fi 
done < ${UNIVERSE} >> ${OUTFILE}

#Clean up
rm -rf ${TMPDIR}
