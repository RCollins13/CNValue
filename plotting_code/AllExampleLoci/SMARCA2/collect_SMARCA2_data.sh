#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required to plot SMARCA2 example locus

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/ ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/
fi
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2 ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2
fi

#PLOT REGION: chr9:1,260,000-2,900,000
chr=9
start=1260000
end=2900000

#Critical regions
distal_start=1500588
distal_end=1713727
proximal_start=2650652
proximal_end=2794964

#####Check for any BRAIN elements in plot & critical regions
cd ${WRKDIR}/data/master_annotations/noncoding/
#Plot region
for file in BRAIN_MASTER*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b ${file} | \
  awk -v OFS="\t" '{ print "chr"$0 }'
done
#Distal critical region
for file in BRAIN_MASTER*; do
  echo "${file}"
  bedtools intersect -wb -a <( echo -e "${chr}\t${distal_start}\t${distal_end}" ) -b ${file} | \
  awk -v OFS="\t" '{ print "chr"$4, $5, $6 }'
done

#####Gather CNV density data for CTRL/GERM CNVs in 5kb bins for DEL/DUP + coding/noncoding
#Write 5kb bins
paste <( seq ${start} 5000 $((${end}-5000)) ) <( seq $((${start}+5000)) 5000 ${end} ) | \
awk -v OFS="\t" '{ print "chr9", $1, $2, "chr9:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.hg19.5kb_bins.bed
#Gather CNVs
for CNV in DEL DUP; do
  for filt in all coding haplosufficient; do
    echo -e "#chr\tstart\tend\tCTRL\tGERM" > \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.${CNV}_${filt}_density.5kb_bins.bed
    sed 's/^chr//g' ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.hg19.5kb_bins.bed | \
    awk -v OFS="\t" '{ print $1, $2, $3 }' | bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E2.GRCh37.${filt}.bed.gz | \
    sed 's/^/chr/g' | awk -v OFS="\t" '{ print $1, $2, $3, $4/38628, $5/63629 }' >> \
    ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.${CNV}_${filt}_density.5kb_bins.bed
  done
done

#####Gather fetal brain RNAseq expression data in 1kb bins for plus/minus strands
#Convert RNAseq data to bedGraphs
for strand in plus minus; do
  ${WRKDIR}/bin/utils/bigWigToBedGraph \
  -chrom=chr7 -start=${start} -end=${end} \
  ${WRKDIR}/data/misc/dataForLocusPlots/fetal_brain_RNAseq_${strand}Strand.bigWig \
  ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_expression_${strand}Strand.bedGraph
done
#Write 1kb bins
paste <( seq ${start} 1000 80359000 ) <( seq 79631000 1000 ${end} ) | \
awk -v OFS="\t" '{ print "chr7", $1, $2, "chr3:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.hg19.1kb_bins.bed
#Gather expression
echo -e "#chr\tstart\tend\tExpression" > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_expression.1kb_bins.bed
awk -v OFS="\t" '{ print $1, $2, $3 }' \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.hg19.1kb_bins.bed | \
bedtools map -o sum -c 4 -a - \
-b ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_expression_plusStrand.bedGraph | \
bedtools map -o sum -c 4 -a - \
-b ${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_expression_minusStrand.bedGraph | \
awk -v OFS="\t" '{ if ($4==".") $4=0; print $0 }' | awk -v OFS="\t" '{ if ($5==".") $5=0; print $0 }' | \
awk -v OFS="\t" '{ print $1, $2, $3, $4+$5 }' >> \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_expression.1kb_bins.bed

#####Gather annotated fetal brain enhancers
fgrep ">" ${WRKDIR}/data/misc/EnhancerAtlas/Fetal_brain.fasta | \
sed -e 's/>//g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/_/\t/g' -e 's/^chr//g' | \
awk -v OFS="\t" '{ if ($4=="") $4="."; print }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -c 4 -o distinct -i - | \
bedtools intersect -wb -a - -b <( echo -e "7\t${start}\t${end}" ) | \
awk -v OFS="\t" '{ print "chr"$1, $2, $3 }' | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.enhancer_elements.bed

#####Gather fetal brain DHS as bedgraph
${WRKDIR}/bin/utils/bigWigToBedGraph \
-chrom=chr7 -start=${start} -end=${end} \
${WRKDIR}/data/misc/dataForLocusPlots/fetal_brain_DHS.bigWig \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.fetal_brain_DHS.bedGraph

#####Gather adult brain H3K27ac as bedgraph
bedtools intersect -wa \
-a adult_brain_H3K27ac_merged.bedGraph -b <( echo -e "chr7\t${start}\t${end}" ) > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.adult_brain_H3K27ac.bedGraph

#####Gather fetal brain H3K4me1 as bedgraph
bedtools intersect -a ${WRKDIR}/data/misc/dataForLocusPlots/fetal_brain_H3K4me1.tagAlign.gz \
-b <( echo -e "chr7\t${start}\t${end}" ) | bedtools genomecov -bga -i - \
-g <( sed 's/^/chr/g' /data/talkowski/rlc47/src/GRCh37.genome ) > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/SMARCA2/fetal_brain_H3K4me1.bedGraph




#Get table of CNVs per phenotype overlapping critical region
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL; do
      bedtools intersect -wb -a <( echo -e "${chr}\t${distal_start}\t${distal_start}" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s | awk -v OFS="\t" '{ print $1, $2, $3, $3/$2 }'
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )
