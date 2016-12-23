#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Ryan's working code for analyses & processing conducted on ERISone (PHC HPC cluster)

#####Set parameters
WRKDIR=/data/talkowski/Samples/rCNVmap
cd ${WRKDIR}
h37=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.fa
TMPDIR=/scratch/miket/rlc47temp/tmp.files/
module load bedtools/2.22.1
module load samtools/1.3
module load hdf/1.8.14
module load anaconda/4.0.5
module rm gcc-4.4
module load gcc/4.9.0
module load bcftools/1.3.1
PICARD=/data/talkowski/tools/bin/picard-tools-1.137/picard.jar
SFARI_ANNO=/data/talkowski/Samples/SFARI/ASC_analysis/annotations

#####Create directory tree
mkdir ${WRKDIR}/lists
mkdir ${WRKDIR}/bin
mkdir ${WRKDIR}/data
mkdir ${WRKDIR}/data/CNV
mkdir ${WRKDIR}/data/CNV/CNV_MASTER
mkdir ${WRKDIR}/data/plot_data
mkdir ${WRKDIR}/data/annotations
mkdir ${WRKDIR}/data/annotations/exonExclusion
mkdir ${WRKDIR}/analysis
mkdir ${WRKDIR}/analysis/BIN_CNV_pileups
mkdir ${WRKDIR}/analysis/BIN_CNV_burdens
mkdir ${WRKDIR}/analysis/BIN_CNV_permutation
mkdir ${WRKDIR}/analysis/Final_Loci
mkdir ${WRKDIR}/analysis/Functional_Enrichments
mkdir ${WRKDIR}/analysis/Functional_Enrichments_exonExclusion

#####Copy CNVs from TAD intolerance Project
#Note: see old code for parameters & syntax used during callset curation
cp /data/talkowski/rlc47/TAD_intolerance/data/CNV/CNV_MASTER/* \
${WRKDIR}/data/CNV/CNV_MASTER/

#####Run TBRden pileup for all tissue types and all phenotypes (coding & noncoding CNVs)
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 100000 -s 25000 \
    -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed \
    -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    /data/talkowski/rlc47/src/GRCh37.genome"
  done
done

#####Run TBRden analysis
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  for CNV in DEL DUP CNV; do
    #Parallelize analyses (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_all"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_noncoding_TBRden_analysis \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.noncoding.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/ \
    ${group}_vs_CTRL_${CNV}_noncoding"
  done
done

#####Run 10k CNV shift permutation tests for all comparisons
#Note: initial p-value cutoff used: 0.05/26802 = 1.865532e-06
#This corresponds to the number of non-overlapping autosomal 100kb bins we tested (after blacklisting N-mask, etc)
for group in DD SCZ DD_SCZ CNCR; do
  echo ${group}
  if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL
  mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split
  for CNV in DEL DUP CNV; do
    echo ${CNV}
    #Coding + noncoding CNVs
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
    #Noncoding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/26802)) print $0 }' > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
    #Split into partitions of 1k permutations each (x100 per comparison)
    for i in $( seq -w 001 100 ); do
      echo ${i}
      #Coding + noncoding CNVs
      if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i} ]; then
        mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}
      fi
      bsub -q normal -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 100000 -N 1000 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_all \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
      #Noncoding CNVs only
      bsub -q normal -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 100000 -N 1000 \
      -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
      -p ${group}_vs_CTRL_${CNV}_noncoding \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.noncoding.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
    done
  done
done

#####Collect results from permutation tests
#Collect results across all permutations per group
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all noncoding; do
      for i in $( seq -w 001 100 ); do
        zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_${filt}.permuted.${i}.bed.gz | \
        fgrep -v "#" | cut -f7 | paste -s
      done > ${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt
      Rscript -e "options(scipen=100);\
      d <- apply(read.table(\"${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt\",header=F),2,sum);\
      write.table(data.frame(\"perms_less_sig\"=100000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/100001),\
        \"${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt\",row.names=F,col.names=T,sep=\"\t\",quote=F)"
      paste <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/001/${group}_vs_CTRL_${CNV}_${filt}.permuted.001.bed.gz | cut -f1-5 ) \
      ${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt > \
      ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed
    done
  done
done
#Compose table of all bins for each comparison
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all noncoding; do
      cat <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -wb -f 1 -r -a - \
      -b ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed | \
      cut -f1-17,23-25 ) \
      <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -v -f 1 -r -a - \
      -b ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed | \
      awk -v OFS="\t" '{ print $0, "NA\tNA\tNA" }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
      cat <( paste <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
        head -n1 | awk '{ print "#"$0 }' ) <( echo -e "perms_less_sig\tperms_as_or_more_sig\tperm_p" ) ) - > \
      ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
    done
  done
done
#Make master table of all bins with summary stats for each bin & each comparison
for dummy in 1; do
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all noncoding; do
        echo -e "${group}.${CNV}.${filt}.obs_p\t${group}.${CNV}.${filt}.perm_p"
      done
    done
  done | paste -s
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all noncoding; do
        zcat ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed | \
        fgrep -v "#" | cut -f17,20
      done | paste - -
    done | paste - - -
  done | paste - - - -
done | paste <( zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed | \
  cut -f1-4 ) - > \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
${WRKDIR}/bin/rCNVmap/bin/summarize_master_burden_file.R \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2
mv ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2 \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
gzip ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
#Print table
# for group in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all noncoding; do
#       echo -e "${group}\n${CNV}\n${filt}"
#       #Count total number of bins tested
#       zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#       sed '1d' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#       #Count number of nominally significant bins from burden test
#       zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#       sed '1d' | awk '{ if ($NF<=0.05) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#       #Count number of Bonferroni significant bins from burden test
#       zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#       sed '1d' | awk '{ if ($NF<=0.05/8875) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#       if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.bed.gz ]; then
#         #Write window results to file
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.bed.gz | \
#         head -n1 | cut -f1-5,8 > ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_windows.bed
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.bed.gz | \
#         sed '1d' | awk '{ if ($NF<=0.05) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | cut -f1-5,8 >> \
#         ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_windows.bed
#         #Write loci results to file
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.bed.gz | \
#         head -n1 | cut -f1-5,8 > ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.bed.gz | \
#         sed '1d' | awk '{ if ($NF<=0.05) print $0 }' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | cut -f1-5,8 | sort -Vk1,1 -k2,2n -k3,3n | \
#         bedtools merge -c 4,5,6 -o distinct,min,min -i - >> ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed
#         #Count number of bins that passed permutation test
#         cat ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_windows.bed | \
#         awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#         #Count number of merged loci that passed permutation test
#         cat ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed | \
#         awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#         #Count raw overlaps with syndromic loci
#         bedtools intersect -f 0.1 -wa -u \
#         -a ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed \
#         -b /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/rCNVs/talkowski_highQual_rCNVs.bed | wc -l
#       else
#         echo -e "TBD\nTBD\nTBD"
#       fi
#     done | paste - - - - - - - - -
#   done
# done
# #Prep for comparison of final loci within and between sets
# for group in DD SCZ DD_SCZ CNCR; do
#   #Make list of all loci
#   for CNV in CNV DEL DUP; do
#     for filt in all noncoding; do
#       fgrep -v "#" ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed
#     done
#   done | sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v group=${group} '{ print $0, group }' > ${TMPDIR}/${group}_master.bed
#   #Make lists of coding/noncoding loci for dels + dups
#   for filt in all noncoding; do
#     for CNV in DEL DUP; do
#       fgrep -v "#" ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed
#     done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/${group}_${filt}.bed
#   done
# done
# #Cat all master lists
# for group in DD SCZ DD_SCZ CNCR; do
#   cat ${TMPDIR}/${group}_master.bed
# done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/all_master.bed
# #Run comparison of final loci between del / dup within a given set
# for group in DD SCZ DD_SCZ CNCR; do
#   echo -e "-\n-"
#   for CNV in DEL DUP; do
#     for filt in all noncoding; do
#       bedtools intersect -r -f 0.5 -c \
#       -a ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed \
#       -b ${TMPDIR}/${group}_${filt}.bed | awk '{ if ($NF==1) print $0 }' | wc -l
#     done
#   done
# done
# #Run comparison of loci specific to phenotype
# for group in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all noncoding; do
#       awk -v group=${group} '{ if ($NF!=group) print $0 }' ${TMPDIR}/all_master.bed | \
#       bedtools intersect -v -r -f 0.5 -b - \
#       -a ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${filt}.${CNV}.significant_loci.bed | wc -l
#     done
#   done
# done
