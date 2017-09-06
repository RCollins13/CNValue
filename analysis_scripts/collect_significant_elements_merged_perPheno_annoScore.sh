#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Collects significant elements per phenotype

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reads arguments
pheno=$1
CNV=$2
VF=$3
filt=$4
list=$5
sig=$6
OUTFILE=$7

#####Prep new tmpdir
TMPDIR=`mktemp -d`

#####Iterates over classes and merges into single master BED file with significant sites
#Allow 50kb distance between significant elements
dist=50000
sigSites_all=`mktemp`
sigSites_master=`mktemp`
while read class path; do
  sigSites=${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.${sig}_sig_elements.bed.gz
  if [ -e ${sigSites} ]; then
    zcat ${sigSites}
  fi
done < ${list} | sort -Vk1,1 -k2,2n -k3,3n | cut -f1-3 > ${sigSites_all}
bedtools merge -i ${sigSites_all} -d ${dist} > ${sigSites_master}
#Rounds master significant sites to floor & ceiling kbs
sigSites_master_round=`mktemp`
Rscript -e "options(scipen=10000); x <- read.table(\"${sigSites_master}\",header=F); \
            x[,2] <- 1000*floor(x[,2]/1000); \
            x[,3] <- 1000*ceiling(x[,3]/1000); \
            write.table(x,\"${sigSites_master_round}\",col.names=F,row.names=F,quote=F,sep=\"\\t\")"

#####Dices sigSites into 1kb intervals and counts number of significant annotations per interval
sigSites_intervals=`mktemp`
while read chr start end; do
  paste <( seq ${start} 1000 $(( ${end}-1000 )) ) <( seq $(( ${start}+1000 )) 1000 ${end} ) | \
  awk -v OFS="\t" -v chr=${chr} '{ print chr, $1, $2 }'
done < ${sigSites_master_round} | sort -Vk1,1 -k2,2n -k3,3n | \
bedtools intersect -c -a - -b ${sigSites_all} > ${sigSites_intervals}

#####Apply Poisson model and return bins with significantly more annotations than expected by chance
sigBins=`mktemp`
Rscript -e "x <- read.table(\"${sigSites_intervals}\",header=F); \
            lambda <- mean(x[,4]); cutoff <- qpois(0.05,lambda,lower.tail=F); \
            write.table(x[which(x[,4]>=cutoff),1:3],\"${sigBins}\",\
            col.names=F,row.names=F,quote=F,sep=\"\\t\")"

#####Retest significant bins and keep those passing a nominal p-value threshold from Fisher's test
sigBins_data=`mktemp`
#Gathers data for test
${WRKDIR}/bin/rCNVmap/bin/gather_annoScore_data.sh -z -N \
-p sigBins \
-o ${sigBins_data} \
${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
${sigBins} \
${h37}
#Get number of cases & controls
nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
         ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
nCTRL=38628
#Runs model
sigBins_retested=`mktemp`
${WRKDIR}/bin/rCNVmap/bin/run_annoScore_model.R \
--onlyfisher \
-o ${sigBins_retested} \
${sigBins_data}.gz \
${nCTRL} \
${nCASE}
#Collects bins passing a nominal threshold from a Fisher's exact test
sigBins_retested_filtered=`mktemp`
awk -v OFS="\t" '{ if ($30<0.05) print $1, $2, $3 }' ${sigBins_retested} > \
${sigBins_retested_filtered}

#####Retest significant bins for all common CNVs and keep those passing a nominal p-value threshold from Fisher's test
sigBins_common_data=`mktemp`
#Gathers data for test
${WRKDIR}/bin/rCNVmap/bin/gather_annoScore_data.sh -z -N \
-p sigBins \
-o ${sigBins_common_data} \
${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.all.bed.gz \
${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E2.GRCh37.all.bed.gz \
${sigBins_retested_filtered} \
${h37}
#Runs model from common data
sigBins_common_retested=`mktemp`
${WRKDIR}/bin/rCNVmap/bin/run_annoScore_model.R \
--onlyfisher \
-o ${sigBins_common_retested} \
${sigBins_common_data}.gz \
${nCTRL} \
${nCASE}

#Collects significant loci from both original and common tests & merges w/25kb buffer
echo -e "#chr\tstart\tend" > ${OUTFILE}
awk -v OFS="\t" '{ if ($30<0.05) print $1, $2, $3 }' ${sigBins_common_retested} | \
bedtools merge -d 25000 -i - >> ${OUTFILE}

#####Clean up
rm -rf ${TMPDIR}