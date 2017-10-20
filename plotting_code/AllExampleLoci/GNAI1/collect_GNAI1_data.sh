#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required to plot GNAI1 example locus

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/ ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/
fi
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1 ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1
fi

#PLOT REGION: chr7:79,630,000-80,360,000

#####Check for any BRAIN elements in plot & critical regions
cd ${WRKDIR}/data/master_annotations/noncoding/
#Plot region
for file in BRAIN_MASTER*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "7\t79630000\t80360000" ) -b ${file} | \
  awk -v OFS="\t" '{ print "chr"$0 }'
done
#Critical region
for file in BRAIN_MASTER*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "7\t80050742\t80197059" ) -b ${file} | \
  awk -v OFS="\t" '{ print "chr"$4, $5, $6 }'
done

#####Gather CNV density data for CTRL/NDD CNVs in 5kb bins for DEL/DUP + coding/noncoding
#Write 5kb bins
paste <( seq 79630000 5000 80355000 ) <( seq 79635000 5000 80360000 ) | \
awk -v OFS="\t" '{ print "chr7", $1, $2, "chr3:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.hg19.5kb_bins.bed
#Gather CNVs
for CNV in DEL DUP; do
  for filt in all; do
    echo -e "#chr\tstart\tend\tCTRL\tNDD" > \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.${CNV}_${filt}_density.5kb_bins.bed
    sed 's/^chr//g' ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.hg19.5kb_bins.bed | \
    awk -v OFS="\t" '{ print $1, $2, $3 }' | bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/NDD/NDD.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    sed 's/^/chr/g' | awk -v OFS="\t" '{ print $1, $2, $3, $4, $5 }' >> \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.${CNV}_${filt}_density.5kb_bins.bed
  done
done

#####Gather fetal brain RNAseq expression data in 1kb bins for plus/minus strands
#Convert RNAseq data to bedGraphs
for strand in plus minus; do
  ${WRKDIR}/bin/utils/bigWigToBedGraph \
  -chrom=chr7 -start=79630000 -end=80360000 \
  ${WRKDIR}/data/misc/dataForLocusPlots/fetal_brain_RNAseq_${strand}Strand.bigWig \
  ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_expression_${strand}Strand.bedGraph
done
#Write 1kb bins
paste <( seq 79630000 1000 80359000 ) <( seq 79631000 1000 80360000 ) | \
awk -v OFS="\t" '{ print "chr7", $1, $2, "chr3:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.hg19.1kb_bins.bed
#Gather expression
echo -e "#chr\tstart\tend\tExpression" > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_expression.1kb_bins.bed
awk -v OFS="\t" '{ print $1, $2, $3 }' \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.hg19.1kb_bins.bed | \
bedtools map -o sum -c 4 -a - \
-b ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_expression_plusStrand.bedGraph | \
bedtools map -o sum -c 4 -a - \
-b ${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_expression_minusStrand.bedGraph | \
awk -v OFS="\t" '{ if ($4==".") $4=0; print $0 }' | awk -v OFS="\t" '{ if ($5==".") $5=0; print $0 }' | \
awk -v OFS="\t" '{ print $1, $2, $3, $4+$5 }' >> \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_expression.1kb_bins.bed

#####Gather annotated fetal brain enhancers
fgrep ">" ${WRKDIR}/data/misc/EnhancerAtlas/Fetal_brain.fasta | \
sed -e 's/>//g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/_/\t/g' -e 's/^chr//g' | \
awk -v OFS="\t" '{ if ($4=="") $4="."; print }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - | \
bedtools intersect -wb -a - -b <( echo -e "7\t79630000\t80360000" ) | \
awk -v OFS="\t" '{ print "chr"$1, $2, $3 }' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.enhancer_elements.bed

#####Gather fetal brain DHS as bedgraph
${WRKDIR}/bin/utils/bigWigToBedGraph \
-chrom=chr7 -start=79630000 -end=80360000 \
${WRKDIR}/data/misc/dataForLocusPlots/fetal_brain_DHS.bigWig \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.fetal_brain_DHS.bedGraph

#####Gather adult brain H3K27ac as bedgraph
bedtools intersect -wa \
-a adult_brain_H3K27ac_merged.bedGraph -b <( echo -e "chr7\t79630000\t80360000" ) > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/GNAI1/GNAI1_locus.adult_brain_H3K27ac.bedGraph




#Get table of CNVs per phenotype overlapping critical region
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL; do
      bedtools intersect -wb -a <( echo -e "7\t80048636\t80196818" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s | awk -v OFS="\t" '{ print $1, $2, $3, $3/$2 }'
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )
