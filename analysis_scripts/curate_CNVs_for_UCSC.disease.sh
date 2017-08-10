#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Converts master CNV bed file to UCSC-format track

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Read arguments
pheno=$1
CNV=$2
VF=$3
context=$4

#####Run
#Write header
echo -e "track name=\"${pheno}_${CNV}_${VF}_${filt}\" itemRgb=\"On\"" > \
${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
#Curate everything except phenotypes
zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
fgrep -v "#" | sed -e 's/DEL/255,0,0/g' -e 's/DUP/0,0,255/g' | awk -v OFS="\t" \
'{ print "chr"$1, $2, $3, "PHENOS", 1000, "+", $2, $3, $5 }' > \
${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt
#Curate phenotypes
zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
fgrep -v "#" | cut -f6 | sed -f ${WRKDIR}/data/misc/HPO_to_phenotype.sed > \
${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.phenos_verbose.txt
Rscript -e "x <- read.table(\"${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.phenos_verbose.txt\")[,1];\
            x <- sapply(x,function(phenos){return(paste(sort(unique(unlist(strsplit(unlist(strsplit(as.character(phenos),split=\",\")),split=\";\")))),collapse=\";\"))});\
            write.table(x,\"${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partB.txt\",col.names=F,row.names=F,quote=F)"
#Combine parts A & B
paste <( cut -f 1-3 ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt ) \
${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partB.txt \
<( cut -f5- ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt ) >> \
${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
#Gzip output
gzip ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt