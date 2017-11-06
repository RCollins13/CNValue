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
#Create master mask of N-masked regions and 100kb flanking telomeres/1Mb flanking centromeres
##############################################################################
cat <( grep -e 'telomere' /data/talkowski/rlc47/src/GRCh37_heterochromatin.bed | \
awk -v OFS="\t" '{ print $1, $2-100000, $3+100000 }' | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' ) \
<( grep -e 'centromere' /data/talkowski/rlc47/src/GRCh37_heterochromatin.bed | \
awk -v OFS="\t" '{ print $1, $2-1000000, $3+1000000 }' | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' ) \
/data/talkowski/rlc47/src/GRCh37_Nmask.bed \
<( grep -e 'X\|Y\|M' ${WRKDIR}/lists/rCNVmap_excluded_loci.CNVs.bed | cut -f1-3 ) | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed 

################
#Run CNV pileups
################
#200kb bins, 10kb steps, 2.5Mb window, ±1Mb window smoothing, CNVs with 100% overlap only
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
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 200000 -s 10000 -d 2500000 -R 2 -I 1 \
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
        bsub -q short -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_pileup_analysis -u nobody \
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
#Filters: ≥250kb, 100kb merge distance, ≥4 protein-coding genes, and <50% SD coverage
minSize=250000
minGenes=4
maxSD=0.5
maxDist=250000
#Note: combine loci significant for DEL or DUP, but not CNV (DEL+DUP)
if [ -e ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered ]; then
  rm -rf ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered
fi
mkdir ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered
for VF in E2 E3 E4 N1; do
  for filt in all coding haplosufficient noncoding intergenic; do
    for CNV in DEL DUP; do
      cat ${WRKDIR}/analysis/large_CNV_segments/master_lists/${CNV}_${VF}_${filt}.signif.bed
    done | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - -d ${maxDist} | \
    awk -v OFS="\t" -v minSize=${minSize} '{ if ($3-$2>=minSize) print $1, $2, $3 }' | \
    # bedtools intersect -c -a - \
    # -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
    # awk -v minGenes=${minGenes} -v OFS="\t" '{ if ($4>=minGenes) print $1, $2, $3 }' | \
    # bedtools coverage -b - \
    # -a ${WRKDIR}/data/master_annotations/noncoding/SegDups.elements.bed | \
    # awk -v maxSD=${maxSD} -v OFS="\t" '{ if ($NF<maxSD) print $1, $2, $3 }' | \
    sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed
  done
done
#Get median size of merged loci loci
VF=E2
filt=all
awk '{ print $3-$2 }' ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed | \
sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'

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
      #     echo -e "chr${chr}\n${start}\n${end}"
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
#Get odds ratios
#NOTE: require CNV to cover >25% of syndromic locus by size
for VF in E2 E3 E4 N1; do
  for filt in all coding haplosufficient noncoding intergenic; do
    for CNV in DEL DUP; do
      #CODE (parallelized below):
      # while read chr start end; do
      #   for dummy in 1; do
      #     echo -e "chr${chr}\n${start}\n${end}"
      #     #iterate over phenos
      #     while read pheno nCASE; do
      #       #Get counts of case/control CNVs
      #       caseCNV=$( bedtools intersect -wb -f 0.25 -a <( echo -e "${chr}\t${start}\t${end}" ) \
      #       -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | wc -l )
      #       controlCNV=$( bedtools intersect -wb -f 0.25 -a <( echo -e "${chr}\t${start}\t${end}" ) \
      #       -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz | wc -l )
      #       controlCNV=$( echo -e "${controlCNV}+1" | bc )
      #       caseNoCNV=$( echo "${nCASE}-${caseCNV}" | bc )
      #       controlNoCNV=$( echo "38628+1-${controlCNV}" | bc )
      #       #Calcluate odds ratio
      #       unset R_HOME
      #       Rscript -e "cat(paste((${caseCNV}/(${controlCNV}+1))/(${caseNoCNV}/(${controlNoCNV}+1))),\"\n\",sep=\"\")"
      #     done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      #       fgrep -v CTRL | cut -f1,8 )
      #   done | paste -s
      # done < ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed > \
      # ${WRKDIR}/analysis/large_CNV_segments/assoc_stats/DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR.bed
      bsub -q short -sla miket_sc -J DEL_DUP_union.${VF}_${filt}.signif.filtered.${CNV}_OR -u nobody \
      "${WRKDIR}/bin/rCNVmap/analysis_scripts/get_ORs_sigCNVsegs.sh ${VF} ${filt} ${CNV}"
    done
  done
done

############################################################
#####Count significant loci per germline & cancer phenotypes
############################################################
VF=E2
filt=all
#Get count (any pheno)
cat ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed | wc -l
#Get count significant in any germline pheno
while read pheno; do
  for CNV in DEL DUP; do
    bedtools intersect -wa -u \
    -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
    -b ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      fgrep -v CTRL | awk '{ if ($2=="GERM") print $1 }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${TMPDIR}/GERM.sig.loci.bed
cat ${TMPDIR}/GERM.sig.loci.bed | wc -l
#Get count significant in any cancer pheno
while read pheno; do
  for CNV in DEL DUP; do
    bedtools intersect -wa -u \
    -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
    -b ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      fgrep -v CTRL | awk '{ if ($2=="CNCR") print $1 }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${TMPDIR}/CNCR.sig.loci.bed
cat ${TMPDIR}/CNCR.sig.loci.bed | wc -l

########################################
#####Locus overlaps vs positive controls
########################################
#Get overlap between germline & cancer significant loci
bedtools intersect \
-a ${TMPDIR}/GERM.sig.loci.bed \
-b ${TMPDIR}/CNCR.sig.loci.bed | wc -l
#Get overlap of germline significant loci & known loci
bedtools intersect -u -a ${TMPDIR}/GERM.sig.loci.bed \
-b ${WRKDIR}/data/master_annotations/noncoding/syndromic_CNVs.elements.bed | wc -l
#Get overlap of cancer significant loci & known loci
bedtools intersect -u -a ${TMPDIR}/CNCR.sig.loci.bed \
-b ${WRKDIR}/data/master_annotations/noncoding/cancer_CNVs.elements.bed | wc -l
#Overlap each pan CNCR-significant locus w/COSMIC genes
VF=E2
filt=all
for CNV in DEL DUP; do
  while read chr start end; do
    for dummy in 1; do
      echo -e "${chr}\t${start}\t${end}"
      TS=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
        -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
        awk '{ print $NF }' | fgrep -wf ${WRKDIR}/data/master_annotations/genelists/COSMIC_tumor_suppressor.genes.list | \
        paste -s -d, )
      if [ -z ${TS} ]; then
        echo "."
      else
        echo "${TS}"
      fi
      ONCO=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
        -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
        awk '{ print $NF }' | fgrep -wf ${WRKDIR}/data/master_annotations/genelists/COSMIC_oncogene.genes.list | \
        paste -s -d, )
      if [ -z ${ONCO} ]; then
        echo "."
      else
        echo "${ONCO}"
      fi
    done | paste -s
  done < <( bedtools intersect -wa -u \
    -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
    -b ${WRKDIR}/analysis/large_CNV_segments/CNCR/CNCR_${CNV}_${VF}_${filt}.signif.bed )
done

##################################################
#####Overlap of genes hit versus constrained genes
##################################################
#Get median count of protein-coding genes per significant locus
VF=E2
filt=all
bedtools intersect -c \
-a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
-b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
cut -f4 | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#Get lists of significant loci per phenotype/CNV combo
for CNV in DEL DUP; do
  #All phenos
  while read pheno; do
    bedtools intersect -wa -u \
    -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
    -b ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
            fgrep -v CTRL | cut -f1 ) | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${TMPDIR}/ALL_PHENOS.${CNV}.sig.loci.bed
  #CNCR or GERM only
  for group in CNCR GERM; do
    while read pheno; do
      bedtools intersect -wa -u \
      -a ${WRKDIR}/analysis/large_CNV_segments/master_lists/filtered/DEL_DUP_union.${VF}_${filt}.signif.filtered.bed \
      -b ${WRKDIR}/analysis/large_CNV_segments/${pheno}/${pheno}_${CNV}_${VF}_${filt}.signif.bed
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v CTRL | awk -v group=${group} '{ if ($2==group) print $1 }' ) | \
              sort -Vk1,1 -k2,2n -k3,3n | uniq > \
    ${TMPDIR}/${group}.${CNV}.sig.loci.bed
  done
done
#Get counts of constrained genes per phenotype/CNV combo
if ! [ -e ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments ]; then
  mkdir ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments
fi
for CNV in DEL DUP; do
  for group in ALL_PHENOS GERM CNCR; do
    paste <( fgrep -wf ${WRKDIR}/data/master_annotations/genelists/ExAC_constrained.genes.list \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed |
      bedtools intersect -c -b - -a ${TMPDIR}/${group}.${CNV}.sig.loci.bed ) \
      <( bedtools intersect -c -a ${TMPDIR}/${group}.${CNV}.sig.loci.bed \
      -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
    awk -v OFS="\t" '{ print $4, $8 }' > \
    ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments/${group}.${CNV}.constrained_genes_count.txt
  done
done 
#Get counts of missense-constrained genes per phenotype/CNV combo
for CNV in DEL DUP; do
  for group in ALL_PHENOS GERM CNCR; do
    paste <( fgrep -wf ${WRKDIR}/data/master_annotations/genelists/ExAC_missense_constrained.genes.list \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed |
      bedtools intersect -c -b - -a ${TMPDIR}/${group}.${CNV}.sig.loci.bed ) \
      <( bedtools intersect -c -a ${TMPDIR}/${group}.${CNV}.sig.loci.bed \
      -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed ) | \
    awk -v OFS="\t" '{ print $4, $8 }' > \
    ${WRKDIR}/analysis/large_CNV_segments/constrained_enrichments/${group}.${CNV}.missense_constrained_genes_count.txt
  done
done 

####################################################################################
#####Get count of cases carrying at least one DEL or DUP at a significant GERM locus
####################################################################################
VF=E2
filt=all
group=GERM
if ! [ -e ${WRKDIR}/${WRKDIR}/analysis/large_CNV_segments/prevalance ]; then
  mkdir ${WRKDIR}/analysis/large_CNV_segments/prevalance
fi
#Iterate over phenotypes & count CNVs (≥25% overlap)
for CNV in DEL DUP; do
  while read pheno nSamp; do
    for dummy in 1; do
      echo "${pheno}"
      bedtools intersect -f 0.25 -wb -a ${TMPDIR}/${group}.${CNV}.sig.loci.bed \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz | \
      cut -f7 | sort | uniq | wc -l
      echo "${nSamp}"
    done | paste -s
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
            sed -n '1,23p' | awk -v OFS="\t" '{ print $1, $NF }' ) > \
  ${WRKDIR}/analysis/large_CNV_segments/prevalance/GERM_${CNV}.prevalance.txt
done















