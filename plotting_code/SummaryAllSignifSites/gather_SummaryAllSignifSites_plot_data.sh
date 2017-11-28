#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to gather all data required for summary of all significant sites

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/plot_data/SummaryAllSignifSites ]; then
  rm -rf ${WRKDIR}/data/plot_data/SummaryAllSignifSites
fi
mkdir ${WRKDIR}/data/plot_data/SummaryAllSignifSites

####Panel a: counts of significant sites by phenotype and DEL/DUP
#Totals of all phenotypes
for wrapper in 1; do
  echo "Any Pheno"
  for CNV in DEL DUP; do
    ##Significant large segments
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | wc -l

    ##Significant genes
    #Count significant genes
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | wc -l

    ##Significant regulatory blocks
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | wc -l
  done
done | paste -s > ${WRKDIR}/data/plot_data/SummaryAllSignifSites/signifSites.countPerPheno.txt
for pheno in GERM NEURO NDD PSYCH NONN; do
  for wrapper in 1; do
    echo -e "${pheno}"
    for CNV in DEL DUP; do
      ##Significant large segments
      #Find column header
      idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | \
             sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
             fgrep "${pheno}_" | fgrep "_p" | cut -f2 )
      #Count significant segments
      sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | \
      awk -v FS="\t" -v idx=${idx} '{ if ($idx<0.05/12480) print NR }' | wc -l

      ##Significant genes
      #Find column header
      idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | \
             sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
             fgrep "${pheno}_" | fgrep "_q" | cut -f2 )
      #Count significant genes
      sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | \
      awk -v FS="\t" -v idx=${idx} '{ if ($idx<0.05) print NR }' | wc -l

      ##Significant regulatory blocks
      #Find column header
      idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
             sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
             fgrep "${pheno}_" | fgrep "_p" | cut -f2 )
      #Count significant genes
      sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
      awk -v FS="\t" -v idx=${idx} '{ if ($idx<0.05) print NR }' | wc -l
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/SummaryAllSignifSites/signifSites.countPerPheno.txt


#####Panel b: Overlap of all sites between phenotypes
#Gather list of loci (either DEL and/or DUP) per pheno
for pheno in GERM NEURO NDD PSYCH NONN; do
  for file in segments genes blocks; do
    if [ -e ${TMPDIR}/all_${pheno}_${file}.txt ]; then
      rm ${TMPDIR}/all_${pheno}_${file}.txt
    fi
  done
  for CNV in DEL DUP; do
    ##Significant large segments
    #Find column header
    idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | \
           sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
           fgrep "${pheno}_" | fgrep "_p" | cut -f2 )
    #Gather significant segments
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_1_2_${CNV}.txt | \
    awk -v FS="\t" -v OFS="\t" -v idx=${idx} '{ if ($idx<0.05/12480) print $1"_"$2"_"$3 }' >> \
    ${TMPDIR}/all_${pheno}_segments.txt
    sort ${TMPDIR}/all_${pheno}_segments.txt | \
    uniq > ${TMPDIR}/all_${pheno}_segments.txt2
    mv ${TMPDIR}/all_${pheno}_segments.txt2 ${TMPDIR}/all_${pheno}_segments.txt

    ##Significant genes
    #Find column header
    idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | \
           sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
           fgrep "${pheno}_" | fgrep "_q" | cut -f2 )
    #Count significant genes
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_3_4_${CNV}.txt | \
    awk -v FS="\t" -v idx=${idx} '{ if ($idx<0.05) print $1 }' >> \
    ${TMPDIR}/all_${pheno}_genes.txt
    sort ${TMPDIR}/all_${pheno}_genes.txt | uniq > ${TMPDIR}/all_${pheno}_genes.txt2
    mv ${TMPDIR}/all_${pheno}_genes.txt2 ${TMPDIR}/all_${pheno}_genes.txt

    ##Significant regulatory blocks
    #Find column header
    idx=$( head -n1 ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
           sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
           fgrep "${pheno}_" | fgrep "_p" | cut -f2 )
    #Count significant genes
    sed '1d' ${WRKDIR}/data/plot_data/suppTables/suppTables_5_6_${CNV}.txt | \
    awk -v FS="\t" -v idx=${idx} '{ if ($idx<0.05) print $4 }' >> \
    ${TMPDIR}/all_${pheno}_blocks.txt
    sort ${TMPDIR}/all_${pheno}_blocks.txt > ${TMPDIR}/all_${pheno}_blocks.txt2
    mv ${TMPDIR}/all_${pheno}_blocks.txt2 ${TMPDIR}/all_${pheno}_blocks.txt
  done
done
#Combine loci across classes
for pheno in GERM NEURO NDD PSYCH NONN; do
  for class in segments genes blocks; do
    cat ${TMPDIR}/all_${pheno}_${class}.txt
  done | sort | uniq > ${TMPDIR}/all_${pheno}_loci.txt
done
#Get master list of all loci
for pheno in GERM NEURO NDD PSYCH NONN; do
  for class in segments genes blocks; do
    cat ${TMPDIR}/all_${pheno}_${class}.txt
  done
done | sort | uniq > ${TMPDIR}/all_master_loci.txt
#Iterate over master list of loci and test for pheno associations
while read ID; do
  assoc=""
  for pheno in GERM NEURO NDD PSYCH NONN; do
    hits=$( fgrep -w ${ID} ${TMPDIR}/all_${pheno}_loci.txt | wc -l )
    if [ ${hits} -gt 0 ]; then
      assoc=$( echo ${assoc} | paste -d_ - <( echo "${pheno}" ) )
    fi
  done
  count=$( echo ${assoc} | sed 's/_/\t/g' | awk '{ print NF }' )
  echo -e "${ID}\t${count}\t${assoc}"
done < ${TMPDIR}/all_master_loci.txt | sed 's/\t_/\t/g' | sort -nrk2,2 > \
${WRKDIR}/data/plot_data/SummaryAllSignifSites/signifSites.phenoSharing.txt






