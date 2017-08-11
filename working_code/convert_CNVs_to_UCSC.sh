#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to clean all CNVs to UCSC track format

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/CNV/UCSC_CNV ]; then
  rm -rf ${WRKDIR}/data/CNV/UCSC_CNV
fi
mkdir ${WRKDIR}/data/CNV/UCSC_CNV

#####Prepare sed file for converting HPOs to long-form phenotypes
while read HPO; do
  awk -v HPO=${HPO} '$5 ~ HPO { print $4 }' \
  ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
  sort | uniq | paste -s -d, | awk -v HPO=${HPO} \
  '{ print "s/"HPO"/"$1"/g" }'
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
cut -f5 | sed 's/\;/\n/g' | sort | uniq | fgrep -v CONTROL ) > \
${WRKDIR}/data/misc/HPO_to_phenotype.sed

#####Iterate over disease phenotypes
while read pheno; do
  #Reinitialize directory if exists
  if [ -e ${WRKDIR}/data/CNV/UCSC_CNV/${pheno} ]; then
    rm -rf ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}
  fi
  mkdir ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}
  #Iterate over filter combinations
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        #CODE (parallelized below):
        # #Write header
        # echo -e "track name=\"${pheno}_${CNV}_${VF}_${filt}\" itemRgb=\"On\"" > \
        # ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
        # #Curate everything except phenotypes
        # zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
        # fgrep -v "#" | sed -e 's/DEL/255,0,0/g' -e 's/DUP/0,0,255/g' | awk -v OFS="\t" \
        # '{ print "chr"$1, $2, $3, "PHENOS", 1000, "+", $2, $3, $5 }' > \
        # ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt
        # #Curate phenotypes
        # zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
        # fgrep -v "#" | cut -f6 | sed -f ${WRKDIR}/data/misc/HPO_to_phenotype.sed > \
        # ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.phenos_verbose.txt
        # Rscript -e "x <- read.table(\"${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.phenos_verbose.txt\")[,1];\
        #             x <- sapply(x,function(phenos){return(paste(sort(unique(unlist(strsplit(unlist(strsplit(as.character(phenos),split=\",\")),split=\";\")))),collapse=\";\"))});\
        #             write.table(x,\"${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partB.txt\",col.names=F,row.names=F,quote=F)"
        # #Combine parts A & B
        # paste <( cut -f 1-3 ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt ) \
        # ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partB.txt \
        # <( cut -f5- ${TMPDIR}/${pheno}.${CNV}.${VF}.${filt}.partA.txt ) >> \
        # ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
        # #Gzip output
        # gzip ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
        #PARALLELIZE w/LSF:
        bsub -q short -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_UCSCformat -u nobody \
        "${WRKDIR}/bin/rCNVmap/analysis_scripts/curate_CNVs_for_UCSC.disease.sh \
        ${pheno} ${CNV} ${VF} ${filt}"
      done
    done
  done
done < <( sed '1d' ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          fgrep -v CTRL | cut -f1 )

#####Repeat same process for controls
for pheno in CTRL; do
  #Reinitialize directory if exists
  if [ -e ${WRKDIR}/data/CNV/UCSC_CNV/${pheno} ]; then
    rm -rf ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}
  fi
  mkdir ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}
  #Iterate over filter combinations
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        #Print status
        echo -e "${pheno} ${CNV} ${VF} ${filt}"
        #Write header
        echo -e "track name=\"${pheno}_${CNV}_${VF}_${filt}\" itemRgb=\"On\"" > \
        ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
        #Curate everything, including phenotype
        zcat ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
        fgrep -v "#" | sed -e 's/DEL/255,0,0/g' -e 's/DUP/0,0,255/g' | awk -v OFS="\t" \
        '{ print "chr"$1, $2, $3, "Healthy_control", 1000, "+", $2, $3, $5 }' | \
        bedtools intersect -f 1 -a - \
        -b <( awk -v OFS="\t" '{ print "chr"$1, 1, $2 }' /data/talkowski/rlc47/src/GRCh37.genome ) >> \
        ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
        #Gzip output
        gzip ${WRKDIR}/data/CNV/UCSC_CNV/${pheno}/${pheno}.${CNV}.${VF}.${filt}.UCSC_track.txt
      done
    done
  done
done



