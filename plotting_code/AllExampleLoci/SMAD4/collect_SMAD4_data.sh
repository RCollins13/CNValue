#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required to plot SMAD4 example locus

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/ ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/
fi
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4 ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4
fi

#####Get BRAIN elements between SMAD4 promoter and disrupted TBR
cd ${WRKDIR}/data/master_annotations/noncoding/
for file in BRAIN_MASTER*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "18\t48550490\t49193871" ) -b ${file} | cut -f4-
done
for file in *CompartmentB*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "18\t48550490\t49193871" ) -b ${file} | cut -f4-
done

#chr18:47365865-53636506

#Write 5kb bins
paste <( seq 44365000 5000 53190000 ) <( seq 44370000 5000 53195000 ) | \
awk -v OFS="\t" '{ print "chr18", $1, $2, "chr3:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/SMAD4_locus.hg19.5kb_bins.bed
#Iterate over tissues and compute compartment scores
while read code tissue; do
  echo ${tissue}
  ${WRKDIR}/bin/utils/bigWigAverageOverBed \
  -bedOut=${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/${tissue}.compartments_5kbBins.bed \
  /data/talkowski/Samples/SFARI/ASC_analysis/annotations/TADs/Schmitt2016/Compartment_primary_cohort/${code}.pc.bw \
  ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/SMAD4_locus.hg19.5kb_bins.bed \
  /dev/null
done < <( echo -e "CO\tcortex\nHC\thippocampus\nimr90\tIMR90\nLG\tlung\nPO\tmuscle\nAO\taorta\nBL\tbladder\nGM12878\tLCL" )

#####Gather CNV density data for CTRL/GERM/CNCR in 5kb bins for DEL/DUP + coding/noncoding
for CNV in DEL DUP; do
  for filt in coding noncoding; do
    echo -e "#chr\tstart\tend\tCTRL\tGERM\tCNCR" > \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/SMAD4_locus.${CNV}_${filt}_density.5kb_bins.bed
    sed 's/^chr//g' ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/SMAD4_locus.hg19.5kb_bins.bed | \
    awk -v OFS="\t" '{ print $1, $2, $3 }' | bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CNCR/CNCR.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    sed 's/^/chr/g' | awk -v OFS="\t" '{ print $1, $2, $3, $4/38628, $5/63629, $6/10844 }' >> \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMAD4/SMAD4_locus.${CNV}_${filt}_density.5kb_bins.bed
  done
done







#CNVs by phenotype for SMAD4 TAD boundary chr18:49130554-49212254
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL DUP; do
      bedtools intersect -wb -a <( echo -e "18\t49130554\t49212254" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )
