#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2018 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Get all target genes for nanostring experiment

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Build table
#Note: must have written CNVs of interest to the following BED file:
# ${TMPDIR}/CNVs_for_nanostring.bed
while read chr start end; do
  for wrapper in 1; do
    #Print coordinates
    echo -e "${chr}\t${start}\t${end}"
    #Print count and IDs for all overlapping genes
    genes=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) \
             -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
             awk '{ print $NF }' | sort | uniq | paste -s -d, )
    echo ${genes} | sed 's/,/\n/g' | sed '/^$/d' | wc -l
    echo ${genes}
    #Print nearest non-overlapping genes (both 3' and 5')
    flanking=$( for subwrapper in 1; do
                  #Downstream
                  bedtools closest -D ref -io -iu \
                  -a <( echo -e "${chr}\t${start}\t${end}" ) \
                  -g /data/talkowski/rlc47/src/GRCh37.genome \
                  -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
                  awk '{ print $(NF-1) }'
                  #Upstream
                  bedtools closest -D ref -io -id \
                  -a <( echo -e "${chr}\t${start}\t${end}" ) \
                  -g /data/talkowski/rlc47/src/GRCh37.genome \
                  -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
                  awk '{ print $(NF-1) }'
                done | paste -s -d, )
    echo ${flanking} | sed 's/,/\n/g' | sed '/^$/d' | wc -l
    if [ -z ${flanking} ]; then
      echo "."
    else
      echo ${flanking}
    fi
    #Print all LCL-expressed genes in overlapping TADs
    TAD_LCL=$( bedtools intersect -wb \
                -a <( echo -e "${chr}\t${start}\t${end}" ) \
                -b <( cat <( cut -f1-3 /data/talkowski/rlc47/TAD_intolerance/data/nuc/TAD/GM12878.TAD.bed ) \
                          <( zcat ${SFARI_ANNO}/TADs/Rao2014/GM12878.TADs.bed.gz | cut -f1-3 ) ) | cut -f4-6 | \
                bedtools intersect -wb -a - \
                         -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
                         awk '{ print $NF }' | sort | uniq | \
                fgrep -wvf ${WRKDIR}/data/master_annotations/genelists/GTEx_Cells_EBVtransformed_lymphocytes_lowly_expressed.genes.list | \
                sort | uniq | paste -s -d, )
    echo ${TAD_LCL} | sed 's/,/\n/g' | sed '/^$/d' | wc -l
    if [ -z ${TAD_LCL} ]; then
      echo "."
    else
      echo ${TAD_LCL}
    fi
    #Print all genes in overlapping TADs
    TAD_all=$( bedtools intersect -wb \
                -a <( echo -e "${chr}\t${start}\t${end}" ) \
                -b <( cat <( cut -f1-3 /data/talkowski/rlc47/TAD_intolerance/data/nuc/TAD/GM12878.TAD.bed ) \
                          <( zcat ${SFARI_ANNO}/TADs/Rao2014/GM12878.TADs.bed.gz | cut -f1-3 ) ) | cut -f4-6 | \
                bedtools intersect -wb -a - \
                         -b ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed | \
                         awk '{ print $NF }' | sort | uniq | paste -s -d, )
    echo ${TAD_all} | sed 's/,/\n/g' | sed '/^$/d' | wc -l
    if [ -z ${TAD_all} ]; then
      echo "."
    else
      echo ${TAD_all}
    fi
  done | paste -s
done < ${TMPDIR}/CNVs_for_nanostring.bed > ${TMPDIR}/genes_for_nanostring.bed

#####Count genes
for wrapper in 1; do
  echo -e "Criteria\ngene_count_all\ngene_count_excluding_1Mb_CNVs_and_16p11.2"
  echo -e "Only genes overlapping CNV"
  awk -v OFS="\n" '{ print $5 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  fgrep -v KCTD13 ${TMPDIR}/genes_for_nanostring.bed | \
  awk -v OFS="\n" '{ if ($3-$2<1000000) print $5 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  echo -e "Genes overlapping CNV + flanking genes"
  awk -v OFS="\n" '{ print $5, $7 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  fgrep -v KCTD13 ${TMPDIR}/genes_for_nanostring.bed | \
  awk -v OFS="\n" '{ if ($3-$2<1000000) print $5, $7 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  echo -e "Genes overlapping CNV + flanking genes + LCL-expressed genes in all overlapping TADs"
  awk -v OFS="\n" '{ print $5, $7, $9 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  fgrep -v KCTD13 ${TMPDIR}/genes_for_nanostring.bed | \
  awk -v OFS="\n" '{ if ($3-$2<1000000) print $5, $7, $9 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  echo -e "Genes overlapping CNV + flanking genes + all genes in all overlapping TADs"
  awk -v OFS="\n" '{ print $5, $7, $9, $11 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
  fgrep -v KCTD13 ${TMPDIR}/genes_for_nanostring.bed | \
  awk -v OFS="\n" '{ if ($3-$2<1000000) print $5, $7, $9, $11 }' ${TMPDIR}/genes_for_nanostring.bed | \
  sed 's/,/\n/g' | sort | uniq | wc -l
done | paste - - -


#########################################
#####Get breakdown of phenotypes per gene
#########################################
#Prep plot data directory, if necessary
if ! [ -e ${WRKDIR}/data/plot_data/zfish_HPO_collection ]; then
  mkdir ${WRKDIR}/data/plot_data/zfish_HPO_collection
fi
#Iterate over CNV classes
for CNV in DEL DUP; do
  #Print header
  for wrapper in 1; do
    echo "gene"
    fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
      fgrep -v "CNCR" | fgrep -v "CTRL" | cut -f4
  done | paste -s > ${WRKDIR}/data/plot_data/zfish_HPO_collection/${CNV}_HPO.txt
  #Get number of CNVs per HPO term 
  for wrapper in 1; do
    echo "CNVdb"
    #Iterate over HPO terms
    while read label HPO; do
      zcat ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.all.bed.gz | \
      awk '{ print $(NF-1) }' | sed 's/\;/\t/g' | fgrep -wf \
      <( echo -e "${HPO}" | sed 's/\;/\n/g' ) | wc -l
    done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
              fgrep -v "CNCR" | fgrep -v "CTRL" | cut -f4-5 )
  done | paste -s >> ${WRKDIR}/data/plot_data/zfish_HPO_collection/${CNV}_HPO.txt
  #Iterate over genes of interest and count CNVs by phenotype
  while read gene; do
    awk -v gene=${gene} '{ if ($4==gene) print $0 }' \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed | \
    bedtools intersect -wa -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.all.bed.gz | \
    awk '{ print $(NF-1) }' > ${TMPDIR}/${gene}_HPOs.txt
    for wrapper in 1; do
      echo "${gene}"
      while read HPO; do
        sed 's/\;/\t/g' ${TMPDIR}/${gene}_HPOs.txt | fgrep -wf \
        <( echo -e "${HPO}" | sed 's/\;/\n/g' ) | wc -l
      done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
                fgrep -v "CNCR" | fgrep -v "CTRL" | cut -f5 )
    done | paste -s >> ${WRKDIR}/data/plot_data/zfish_HPO_collection/${CNV}_HPO.txt
  done < <( echo -e "USP7\nADAMTS2\nDMRT2\nATG5\nPCNX\nUNCX\nSNTG1\nTGM6\nLINGO1\n\
                     SEMA5A\nRELL1\nROCK1\nCREBBP\nNRXN1\nCTDP1\nFAM19A5\nNTN1\n\
                     NUP155\nC16orf72\nCOLEC12\nDOK6\nERC2\nMARCH1\nMEI4\nTAC1\n\
                     TPO\nLRPPRC\nMTX2\nSUCLG1\nUGT2B7\nCOL23A1" )
done








