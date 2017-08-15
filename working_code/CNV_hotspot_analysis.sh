#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to identify large CNV hotspots (>500kb) for rare CNV project

###############
#Set parameters
###############
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

##############################################################################
#Create master mask of N-masked regions and 1Mb flanking telomeres/centromeres
##############################################################################
grep -e 'centromere\|telomere' /data/talkowski/rlc47/src/GRCh37_heterochromatin.bed | \
awk -v OFS="\t" '{ print $1, $2-1000000, $3+1000000 }' | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' | \
cat - /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
<( grep -e 'X\|Y\|M' ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed 

################
#Run CNV pileups
################
#250kb bins, 25kb steps
#Clear directory
if [ -e ${WRKDIR}/analysis/BIN_CNV_pileups ]; then
  rm -rf ${WRKDIR}/analysis/BIN_CNV_pileups
fi
mkdir ${WRKDIR}/analysis/BIN_CNV_pileups
#Iterate over phenotypes & CNV combos
while read pheno; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        #Refresh directory
        if [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno} ]; then
          rm -rf ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}
        fi
        mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}
        #Parallelize intersections
        bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_${filt}_binned_pileup \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 250000 -s 25000 -d 1000000 -r 0 \
        -o ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}/${pheno}.${CNV}.${VF}.${filt}.BIN_CNV_pileup.bed \
        -x ${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed  \
        ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
        /data/talkowski/rlc47/src/GRCh37.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 )

###########################################################
#Run statistical analysis of case vs control pileup burdens
###########################################################
#Initialize directory
if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens ]; then
  rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens
fi
mkdir ${WRKDIR}/analysis/BIN_CNV_burdens
#Iterate over all phenotypes
while read pheno color; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}
  #Iterate over all CNV classes
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        #Parallelize analyses (LSF)
        bsub -q normal -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_pileup_analysis -u nobody \
        -o ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${smooth}kb_smoothed/${pheno}_${CNV}_${filter}.out \
        -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${smooth}kb_smoothed/${pheno}_${CNV}_${filter}.err \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
        ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/CTRL.${CNV}.${VF}.${filt}.BIN_CNV_pileup.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}/${pheno}.${CNV}.${VF}.${filt}.BIN_CNV_pileup.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/ \
        ${pheno}_${CNV}_${VF}_${filt} 0.00000001 ${color}"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v CTRL | cut -f1,7 )

#############################################
#Collect & merge significant loci per disease
#############################################
#Set significance threshold
sig=0.00000001
#Initialize directory
if [ -e ${WRKDIR}/analysis/large_CNV_segments ]; then
  rm -rf ${WRKDIR}/analysis/large_CNV_segments
fi
mkdir ${WRKDIR}/analysis/large_CNV_segments
#Iterate over phenotypes
while read pheno; do
  if [ -e ${WRKDIR}/analysis/large_CNV_segments/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/large_CNV_segments/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/large_CNV_segments/${pheno}
  #Iterate over all CNV classes
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_${CNV}_${VF}_${filt}.TBRden_results.bed.gz | \
        awk -v sig=${sig} -v OFS="\t" '{ if ($NF<sig) print $1, $2, $3 }' | \
        bedtools merge -i - > \
        ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v CTRL | cut -f1 )

#####################################################
#Collect & merge significant loci across all diseases
#####################################################
#Get master list of all significant loci
if [ -e ${WRKDIR}/analysis/large_CNV_segments/master_lists ]; then
  rm -rf ${WRKDIR}/analysis/large_CNV_segments/master_lists
fi
mkdir ${WRKDIR}/analysis/large_CNV_segments/master_lists
for CNV in CNV DEL DUP; do
  for VF in E2 E3 E4 N1; do
    for filt in all coding haplosufficient noncoding intergenic; do
      while read pheno; do
        cat ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
      done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      fgrep -v CTRL | cut -f1 ) | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
      ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed
    done
  done
done

###############################
#Filter master significant loci 
###############################
#Filters: ≥500kb, ≥4 protein-coding genes, and <30% SD coverage
minSize=500000
minGenes=4
maxSD=0.3
#Note: combine loci significant for DEL or DUP, but not CNV (DEL+DUP)
if [ -e ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered ]; then
  rm -rf ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered
fi
mkdir ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered
for VF in E2 E3 E4 N1; do
  for filt in all coding haplosufficient noncoding intergenic; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed
    done | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - | \
    awk -v OFS="\t" -v size=${size} '{ if ($3-$2>=size) print $1, $2, $3 }' | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    awk -v minGenes=${minGenes} -v OFS="\t" '{ if ($4>=minGenes) print $1, $2, $3 }' | \
    bedtools coverage -b - \
    -a ${WRKDIR}/data/master_annotations/noncoding/SegDups.elements.bed | \
    awk -v maxSD=${maxSD} -v OFS="\t" '{ if ($NF<maxSD) print $1, $2, $3 }' | \
    sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed
  done
done

########################################################################
#Get association statistics per phenotype per significant filtered locus
########################################################################
if [ -e ${WRKDIR}/analysis/large_CNV_segments/assoc_stats ]; then
  rm -rf ${WRKDIR}/analysis/large_CNV_segments/assoc_stats
fi
mkdir ${WRKDIR}/analysis/large_CNV_segments/assoc_stats
#Get minimum p-values
for VF in E2 E3 E4 N1; do
  for filt in all coding haplosufficient noncoding intergenic; do
    for CNV in DEL DUP; do
      #CODE (parallelized below):
      # while read chr start end; do
      #   for dummy in 1; do
      #     echo -e "chr${chr}n${start}\n${end}"
      #     #iterate over phenos
      #     while read pheno; do
      #       bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
      #       -b <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${pheno}_${CNV}_${VF}_${filt}.TBRden_results.bed.gz | \
      #       sed '1d' ) | awk '{ print $NF }' > ${TMPDIR}/${pheno}_${CNV}_${VF}_${filt}_chr${chr}_${start}_${end}.pvals.tmp
      #       unset R_HOME
      #       Rscript -e "cat(paste(format(min(read.table(\"${TMPDIR}/${pheno}_${CNV}_${VF}_${filt}_chr${chr}_${start}_${end}.pvals.tmp\",header=F)[,1]),scientific=T),\"\n\",sep=\"\"))"
      #     done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      #       fgrep -v CTRL | cut -f1 )
      #   done | paste -s
      # done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed > \
      # ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal.bed 
      #PARALLELIZE:
      bsub -q short -sla miket_sc -J DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_pVal -u nobody \
      "${WRKDIR}/bin/rCNVmap/analysis_scripts/get_lowest_pVals_sigCNVsegs.sh ${VF} ${filt} ${CNV}"
    done
  done
done











