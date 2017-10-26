#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Generates a matrix of jaccard distances based on CNVs overlapping pairs of BED-style loci

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reads arguments
elements=$1   #chr, start, end, ID
CNVs=$2       #chr, start, end, ID
OUTFILE=$3    #Matrix of jaccard indexes

# #DEV PARAMETERS
# elements=/data/talkowski/Samples/rCNVmap/analysis/perAnno_burden/signif_elements/all_merged/all_classes_haplosuffDELnoncodingDUP_E4.signif_loci.merged.filtered.bed
# CNVs=/scratch/miket/rlc47temp/tmp.files/GERM_CNCR_CNVs_pooled.bed
# OUTFILE=/scratch/miket/rlc47temp/tmp.files/test_jaccard_out.txt

#####Make new temporary directory
TMPDIR=`mktemp -d`

#####Write list of CNV-element pairs
bedtools intersect -wa -wb -a ${elements} -b ${CNVs} | \
cut -f4,8 | sort -k1,1 -k2,2 | uniq > \
${TMPDIR}/element_CNV_pairs.txt

#####Write header to outfile
cut -f4 ${elements} | paste -s | awk -v OFS="\t" '{ print "IDs", $0 }' > ${OUTFILE}

#####Iterate over each element & write jaccards to 
while read chrA startA endA IDA; do
  for wrapper in 1; do
    echo ${IDA}
    #####Iterate over each other element
    while read chrB startB endB IDB; do
      #Check if elements are on the same chromosome and within 5Mb
      if [ ${chrA} == ${chrB} ] && [ $(( ${startB} - ${endA} )) -le 5000000 ]; then
        #Check if any CNVs hit both
        both=$( fgrep -w ${IDA} ${TMPDIR}/element_CNV_pairs.txt | cut -f2 | \
                fgrep -wf - ${TMPDIR}/element_CNV_pairs.txt | fgrep -w ${IDB} | wc -l )
        #Report zero if no overlap
        if [ ${both} -eq 0 ]; then
          echo 0
          #Otherwise calculate Jaccard distance
        else
          either=$( grep -e "${IDA}\b\|${IDB}\b" ${TMPDIR}/element_CNV_pairs.txt | \
                    cut -f2 | sort | uniq | wc -l )
          echo -e "scale=8;(( ${both}/${either} ))" | bc
        fi
      else
        echo 0
      fi
    done < ${elements}
  done | paste -s
done < ${elements} >> ${OUTFILE}

#####Clean up
rm -rf ${TMPDIR}