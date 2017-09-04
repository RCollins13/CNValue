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
sigBins_merged=`mktemp`
bedtools merge -i ${sigBins} -d ${dist} > ${sigBins_merged}

#####Retest merged significant bins and keep those passing a nominal p-value threshold
sigBins_merged_data=`mktemp`
#Gathers data for test
${WRKDIR}/bin/rCNVmap/bin/gather_annoScore_data.sh -z \
-p sigBins_merged \
-o ${sigBins_merged_data} \
${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
${sigBins_merged} \
${h37}
#Get number of cases & controls
nCASE=$( awk -v pheno=${pheno} '{ if ($1==pheno) print $4 }' \
         ${WRKDIR}/data/plot_data/figure1/sample_counts_by_group.txt )
nCTRL=38628
#Runs model
${WRKDIR}/bin/rCNVmap/bin/run_annoScore_model.R \
-o ${sigBins_merged_retested} \
${annoScoreData} \
${nCTRL} \
${nCASE}
    #Compress output
    if [ -e ${OUTFILE} ]; then
      gzip -f ${OUTFILE}
    fi
  else
    echo "OUTPUT FOR ${pheno} ${CNV} ${VF} ${filt} vs ${path} IS MISSING. SKIPPING..."
  fi
done < ${list}






#####Iterates over elements list and runs model
while read class path; do
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}
  fi
  if ! [ -e ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class} ]; then
    mkdir ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}
  fi
  STATS=${WRKDIR}/analysis/perAnno_burden/${pheno}/${CNV}/${VF}/${filt}/${pheno}.${CNV}.${VF}.${filt}.${class}.annoScore_stats.bed.gz
  if [ -e ${STATS} ]; then
    #Nominally significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($44<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.nom_sig_elements.bed
    #BH-corrected significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($45<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.BH_sig_elements.bed
    #Bonferroni-corrected significant in cases
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($46<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.bonf_sig_elements.bed
    #Gzip all
    gzip -f ${WRKDIR}/analysis/perAnno_burden/signif_elements/${pheno}/${class}/${pheno}.${CNV}.${VF}.${filt}.${class}.*_sig_elements.bed*

    #Nominally significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($47<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.nom_sig_elements.bed
    #BH-corrected significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($48<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.BH_sig_elements.bed
    #Bonferroni-corrected significant in controls
    zcat ${STATS} | fgrep -v "#" | awk -v OFS="\t" '{ if ($49<0.05) print $1, $2, $3, $4 }' > \
    ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.bonf_sig_elements.bed
    #Gzip all
    gzip -f ${WRKDIR}/analysis/perAnno_burden/signif_elements/CTRL/${class}/${pheno}_vs_CTRL.${CNV}.${VF}.${filt}.${class}.*_sig_elements.bed*
  else
    echo "OUTPUT FOR ${pheno} ${CNV} ${VF} ${filt} vs ${path} IS MISSING. SKIPPING..."
  fi
done < ${list}


#####Clean up
rm ${sigSites_master}