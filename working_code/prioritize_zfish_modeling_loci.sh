#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2018 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Prioritization of rCNV loci for zebrafish modeling

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Create directory if necessary
if ! [ -e ${WRKDIR}/analysis/zfish_modeling ]; then
  mkdir ${WRKDIR}/analysis/zfish_modeling
fi

#####Create master table of criteria for prioritizing gene duplications
#Set criteria
VF=E4; context=exonic; CNV=DUP
#Write header
for wrapper in 1; do
  echo -e "#gene\tCTRL_DUP\tGERM_DUP\tn_DUP_sig_pheno\tDUP_sig_phenos"
  echo -e "best_pheno\tbest_pheno_DUP\tbest_pheno_OR\tbest_pheno_Zscore\tbest_pheno_FDRq"
  echo -e "n_DEL_sig_pheno\tDEL_sig_phenos\tCNCR_DEL_sig\tCNCR_DUP_sig"
  echo -e "Brain_expressed\tn_organs_highly_expressed\tn_organs_lowly_expressed\tn_organs_not_expressed"
  echo -e "ExAC_LoF_constrained\tExAC_LoF_tolerant\tExAC_Mis_constrained\tExAC_Mis_tolerant\tRVIS_pct"
  echo -e "DDG2P_dominant_LoF\tDDG2P_dominant_other"
  echo -e "n_zfish_orthologues\tzfish_orthologue_IDs"
  echo -e "n_SSC_proband_DUP\tn_SSC_sibling_DUP"
  echo -e "extTADA_disease_assoc\tAutosomal_dominant_disease\tClinVar_disease_assoc\tOMIM_disease_assoc\tClinGen_HI\tSFARIGene\tSanders_TADA"
  echo -e "Coe2017_LGD\tCoe2017_MIS"
done | paste -s > ${WRKDIR}/analysis/zfish_modeling/DUP_sig_genes.wAnnotations.txt
#Iterate over genes and collect all criteria
while read gene; do
  for wrapper in 1; do

    #Print gene ID
    echo "${gene}"

    #Write pheno association stats to TMP file
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      awk -v OFS="\t" -v gene=${gene} -v pheno=${pheno} \
      '{ if ($1==gene) print pheno, $8, $14, $26, $41, $44 }' \
      ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_${context}.geneScore_stats.txt
    done > ${TMPDIR}/${gene}.pheno_stats.txt

    #Print total number of DUP in control & case
    fgrep -w GERM ${TMPDIR}/${gene}.pheno_stats.txt | cut -f2,3

    #Print number of significant phenos & list them
    awk '{ if ($NF<0.05) print $1 }' ${TMPDIR}/${gene}.pheno_stats.txt | wc -l
    awk '{ if ($NF<0.05) print $1 }' ${TMPDIR}/${gene}.pheno_stats.txt | paste -s -d,

    #Print association stats for strongest p-value
    sort -nk6,6 ${TMPDIR}/${gene}.pheno_stats.txt | head -n1 | cut --complement -f2

    #Count number of DEL-sig phenos & list them
    for pheno in GERM NEURO NDD PSYCH SOMA; do
      awk -v OFS="\t" -v gene=${gene} -v pheno=${pheno} \
      '{ if ($1==gene) print pheno, $8, $14, $26, $41, $44 }' \
      ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_DEL_${VF}_${context}.geneScore_stats.txt
    done > ${TMPDIR}/${gene}.pheno_stats.DEL.txt
    awk '{ if ($NF<0.05) print $1 }' ${TMPDIR}/${gene}.pheno_stats.DEL.txt | wc -l
    DEL_sig=$( awk '{ if ($NF<0.05) print $1 }' ${TMPDIR}/${gene}.pheno_stats.DEL.txt | paste -s -d, )
    if [ -z ${DEL_sig} ]; then
      DEL_sig="NA"
    fi
    echo "${DEL_sig}"

    #DEL & DUP sig in cancer
    pheno=CNCR
    for CNV in DEL DUP; do
      awk -v OFS="\t" -v gene=${gene} '{ if ($1==gene && $44<0.05) print 1 }' \
      ${WRKDIR}/analysis/perGene_burden/${pheno}/${pheno}_DEL_${VF}_${context}.geneScore_stats.txt | wc -l
    done
    CNV=DUP

    #Expressed in brain
    echo "${gene}" | fgrep -wvf \
    ${WRKDIR}/data/master_annotations/genelists/BRAIN_MASTER_INTERSECTION.Lowly_Expressed.genes.list | wc -l

    #Number of organs highly expressed
    echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
    ${WRKDIR}/data/master_annotations/genelists/*_MASTER_UNION.Highly_Expressed.genes.list | wc -l    

    #Number of organs lowly expressed
    echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
    ${WRKDIR}/data/master_annotations/genelists/*_MASTER_INTERSECTION.Lowly_Expressed.genes.list | wc -l

    #Number of organs not expressed
    echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
    ${WRKDIR}/data/master_annotations/genelists/*_MASTER_INTERSECTION.Not_Expressed.genes.list | wc -l

    #ExAC constraint
    for set in ExAC_constrained ExAC_haplosufficient ExAC_missense_constrained ExAC_missense_tolerant; do
      echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
      ${WRKDIR}/data/master_annotations/genelists/${set}.genes.list | wc -l
    done

    #RVIS percentile
    RVIS=$( echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
    ${WRKDIR}/data/misc/RVIS_Unpublished_ExAC_May2015.txt | awk '{ print $(NF-2) }' )
    if [ -z ${RVIS} ]; then
      RVIS="NA"
    fi
    echo ${RVIS}

    #DDG2P dominant LOF
    echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
    ${WRKDIR}/data/master_annotations/genelists/DDG2P*Dominant_LOF*genes.list | \
    cut -f2 -d\: | sort | uniq | wc -l

    #DD dominant GOF/other
    for mech in GOF Unknown; do
      cat ${WRKDIR}/data/master_annotations/genelists/DDG2P*Dominant_${mech}*genes.list
    done | sort | uniq | fgrep -wf <( echo -e "${gene}" | sed 's/\-/_/g' ) | wc -l

    #Zfish orthologues (count and list)
    orthos=$( awk -v FS="\t" -v gene=${gene} '{ if ($5==gene) print $8 }' \
    ${WRKDIR}/data/misc/mapped_zfish_human_orthologues.dup_sig_genes.txt | \
    sort | uniq | paste -s -d\, )
    if [ -z ${orthos} ]; then
      orthos="NA"
    fi
    echo "${orthos}" | sed 's/,/\n/g' | sort | uniq | wc -l
    echo ${orthos}

    #SFARI proband & sib CNVs (count)
    for mem in probands siblings; do
      fgrep -w $( echo ${gene} | sed 's/\-/_/g' ) \
      <( sed 's/\-/_/g' ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed ) | \
      bedtools intersect -u -b - \
      -a ${TMPDIR}/SSC_CNVs.p10E9.hg19.${mem}.${CNV}.bed | \
      cut -f4 | sort | uniq | wc -l
    done

    #Other disease association lists
    for set in extTADA Autosomal_dominant_disease ClinVar_disease_associated \
    Germline_disease_HPO_associated ClinGen_haploinsufficient_low_confidence \
    ASD_SFARIGene ASD_TADA; do
      hit=$( echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
      ${WRKDIR}/data/master_annotations/genelists/${set}*.genes.list | wc -l )
      if [ ${hit} -gt 1 ]; then
        hit=1
      fi
      echo ${hit}
    done

    #Coe 2017 LGD & MIS
    for mech in LGD MIS; do
      echo -e "${gene}" | sed 's/\-/_/g' | fgrep -wf - \
      ${WRKDIR}/data/master_annotations/genelists/Coe2017_ASD_ID*${mech}.genes.list | wc -l
    done

  done | paste -s
done < ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_FINAL_sig.genes.list >> \
${WRKDIR}/analysis/zfish_modeling/DUP_sig_genes.wAnnotations.txt
