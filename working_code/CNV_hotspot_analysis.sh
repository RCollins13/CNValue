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
while read pheno; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno} ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}
  #Iterate over all CNV classes
  for CNV in DEL DUP CNV; do
    for VF in E2 E3 E4 N1; do
      for filt in all coding haplosufficient noncoding intergenic; do
        #Parallelize analyses (LSF)
        bsub -q normal -sla miket_sc -J ${pheno}_${CNV}_${VF}_${filt}_pileup_analysis -u nobody \
        -o ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${smooth}kb_smoothed/${pheno}_${CNV}_${filter}.out \
        -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/${smooth}kb_smoothed/${pheno}_${CNV}_${filter}.err \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
        ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
        ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}/ \
        ${pheno}_${CNV}_${VF}_${filt} 0.00000001 green"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v CTRL | cut -f1 )

# #####Get significant urCNV loci from 100kb smoothed pileups, allowing for up to 50kb between bins that are merged
# for CNV in CNV DEL DUP; do
#   echo -e "\n\n${CNV}\n\n"
#   while read pheno eti tier descrip include exclude color n; do
#     echo -e "${pheno}\t${n}\t${descrip}"
#     for filter in all coding dispensable noncoding intergenic; do
#       if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/urCNV_100kb_smoothed/${pheno}_vs_CTRL_${CNV}_${filter}.TBRden_results.bed.gz ]; then
#         zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/urCNV_100kb_smoothed/${pheno}_vs_CTRL_${CNV}_${filter}.TBRden_results.bed.gz | \
#         awk '{ if ($NF<=0.0000009545993) print $0 }' | bedtools merge -d 50000 -i - | wc -l
#       else
#         echo "NA"
#       fi
#     done
#   done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
#     fgrep -v "CTRL" ) | paste - - - - - -
# done

# # #####Run TBRden analysis (COE METHOD)
# # #rCNVs
# # while read pheno eti tier descrip include exclude color n; do
# #   if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL ]; then
# #     rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL
# #   fi
# #   mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL
# #   for smooth in 0 50 100; do
# #     if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed ]; then
# #       rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed
# #     fi
# #     mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed
# #     for CNV in DEL DUP CNV; do
# #       for filter in all coding dispensable noncoding intergenic; do
# #         #Parallelize analyses (LSF)
# #         bsub -q normal -sla miket_sc -J ${pheno}_${CNV}_TBRden_analysis -u nobody \
# #         -o ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed/${pheno}_vs_CTRL_${CNV}_${filter}.out \
# #         -e ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed/${pheno}_vs_CTRL_${CNV}_${filter}.err \
# #         "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
# #         ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/${smooth}kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
# #         ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}/${smooth}kb_smoothed/${pheno}.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
# #         ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${smooth}kb_smoothed/ \
# #         ${pheno}_vs_CTRL_${CNV}_${filter} 0.0000003818455 ${color}"
# #       done
# #     done
# #   done
# # done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" )





















# #####Run 1k CNV shift direct permutation tests for all comparisons
# #Note: changed from old 100k matched Fisher permutation
# #Note: initial p-value cutoff used: 0.05/52378 = 9.545993e-07 = 0.0000009545993
# #This corresponds to the number of non-overlapping 50kb bins we tested (after blacklisting N-mask, etc)
# #50kb chosen since minimum CNV size = 50kb, so max # independent tests = size of genome / 50kb
# for pheno in DD SCZ DD_SCZ CNCR; do
#   echo ${pheno}
#   if [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL ]; then
#     rm -rf ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL
#   fi
#   mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL
#   mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split
#   for CNV in DEL DUP CNV; do
#     echo ${CNV}
#     #Coding + noncoding CNVs
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_all.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_all.Bonferroni.bed
#     bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.all.1k_permute -u nobody \
#     "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -N 1000 -t upper \
#     -p ${pheno}_vs_CTRL_${CNV}_all \
#     ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_all.Bonferroni.bed \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/"
#     #Coding CNVs only
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_coding.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_coding.Bonferroni.bed
#     bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.coding.1k_permute -u nobody \
#     "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -c -N 1000 -t upper \
#     -p ${pheno}_vs_CTRL_${CNV}_coding \
#     -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#     ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_coding.Bonferroni.bed \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/"
#     #Noncoding CNVs only
#     zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
#     awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
#     bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.noncoding.1k_permute -u nobody \
#     "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -n -N 1000 -t upper \
#     -p ${pheno}_vs_CTRL_${CNV}_noncoding \
#     -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#     ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed \
#     ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/"
#   done
# done
# # #OLD CODE USED FOR FISHER PERMUTATION
# #Split into partitions of 1k permutations each (x100 per comparison)
# # for i in $( seq -w 001 100 ); do
# #   echo ${i}
# #   #Coding + noncoding CNVs
# #   if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i} ]; then
# #     mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i}
# #   fi
# #   bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
# #   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 1000000 -N 100 \
# #   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i}/${pheno}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
# #   -p ${pheno}_vs_CTRL_${CNV}_all \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_all.Bonferroni.bed"
# #   #Coding CNVs only
# #   bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
# #   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
# #   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i}/${pheno}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
# #   -p ${pheno}_vs_CTRL_${CNV}_coding \
# #   -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
# #   #Noncoding CNVs only
# #   bsub -q short -sla miket_sc -J ${pheno}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
# #   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
# #   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i}/${pheno}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
# #   -p ${pheno}_vs_CTRL_${CNV}_noncoding \
# #   -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}.${CNV}.GRCh37.bed.gz \
# #   ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
# # done


# #####Collect results from permutation tests
# #Collect results across all permutations per pheno
# # for pheno in DD SCZ DD_SCZ CNCR; do
# #   for CNV in CNV DEL DUP; do
# #     for filt in all coding noncoding; do
# #       for i in $( seq -w 001 100 ); do
# #         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/${i}/${pheno}_vs_CTRL_${CNV}_${filt}.permuted.${i}.bed.gz | \
# #         fgrep -v "#" | cut -f7 | paste -s
# #       done > ${TMPDIR}/${pheno}_${CNV}_${filt}_perm.in.txt
# #       Rscript -e "options(scipen=100);\
# #       d <- apply(read.table(\"${TMPDIR}/${pheno}_${CNV}_${filt}_perm.in.txt\",header=F),2,sum);\
# #       write.table(data.frame(\"perms_less_sig\"=10000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/10001),\
# #         \"${TMPDIR}/${pheno}_${CNV}_${filt}_perm.out.txt\",row.names=F,col.names=T,sep=\"\t\",quote=F)"
# #       paste <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/perm_split/001/${pheno}_vs_CTRL_${CNV}_${filt}.permuted.001.bed.gz | cut -f1-5 ) \
# #       ${TMPDIR}/${pheno}_${CNV}_${filt}_perm.out.txt > \
# #       ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed
# #     done
# #   done
# # done
# #Compose table of all bins for each comparison
# for pheno in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       cat <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#       sed '1d' | bedtools intersect -wb -f 1 -r -a - \
#       -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
#       cut -f1-17,26 ) \
#       <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#       sed '1d' | bedtools intersect -v -f 1 -r -a - \
#       -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
#       awk -v OFS="\t" '{ print $0, "NA" }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
#       cat <( paste <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
#         head -n1 | awk '{ print "#"$0 }' ) <( echo -e "perm_sim_p" ) ) - > \
#       ${WRKDIR}/analysis/Final_Loci/${pheno}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
#       gzip -f ${WRKDIR}/analysis/Final_Loci/${pheno}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
#     done
#   done
# done
# #Make master table of all bins with summary stats for each bin & each comparison
# for dummy in 1; do
#   for pheno in CTRL DD SCZ CNCR; do
#     for CNV in DEL DUP; do
#       zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
#       fgrep -v "#" | cut -f5 > ${TMPDIR}/${pheno}.${CNV}.all.CNVcounts.tmp
#       for filt in coding noncoding; do
#         zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
#         fgrep -v "#" | cut -f5 > ${TMPDIR}/${pheno}.${CNV}.${filt}.CNVcounts.tmp
#         echo ${TMPDIR}/${pheno}.${CNV}.${filt}.CNVcounts.tmp
#       done
#     done
#   done
#   for pheno in DD SCZ DD_SCZ CNCR; do
#     for CNV in CNV DEL DUP; do
#       for filt in all coding noncoding; do
#         zcat ${WRKDIR}/analysis/Final_Loci/${pheno}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed.gz | \
#         fgrep -v "#" | cut -f17-18 > ${TMPDIR}/${pheno}.${CNV}.${filt}.pvals.tmp
#         echo ${TMPDIR}/${pheno}.${CNV}.${filt}.pvals.tmp
#       done
#     done
#   done
# done > ${TMPDIR}/list_of_paths.txt
# for dummy in 1; do
#   zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
#   cut -f1-4 | head -n1
#   for pheno in CTRL DD SCZ CNCR; do
#     for CNV in DEL DUP; do
#       for filt in coding noncoding; do
#         echo "${pheno}.${CNV}.${filt}.count"
#       done
#     done
#   done
#   for pheno in DD SCZ DD_SCZ CNCR; do
#     for CNV in CNV DEL DUP; do
#       for filt in all coding noncoding; do
#         echo -e "${pheno}.${CNV}.${filt}.obs_p\n${pheno}.${CNV}.${filt}.perm_p"
#       done
#     done
#   done
# done | paste -s > ${TMPDIR}/header.txt
# paste <( zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
#   fgrep -v "#" | cut -f1-4 ) $( cat ${TMPDIR}/list_of_paths.txt ) | \
# cat ${TMPDIR}/header.txt - > ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
# gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
# ${WRKDIR}/bin/rCNVmap/bin/summarize_master_burden_file.R \
# ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz \
# ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2
# mv ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2 \
# ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
# gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed

# #####Cut master table of bins that are significant per disease 
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   echo ${pheno}
#   if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${pheno} ]; then
#     rm -rf ${WRKDIR}/analysis/Final_Loci/significant/${pheno}
#   fi
#   mkdir ${WRKDIR}/analysis/Final_Loci/significant/${pheno}
#   for CNV in ANY_CNV CNV DEL DUP; do
#     echo ${CNV}
#     for filt in ANY_FILTER all coding noncoding; do
#       echo ${filt}
#       #Parallelize:
#       bsub -q short -sla miket_sc -J ${pheno}_${CNV}_${filt}_filterMasterBurden \
#       "${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file_parallelized.sh \
#       ${pheno} ${CNV} ${filt}"
#       # list=`mktemp`
#       # echo -e "${pheno}.${CNV}.${filt}.obs_p" > ${list}
#       # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
#       # -o ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_nom_signif_bins.bed \
#       # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
#       # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R -t 0.05 \
#       # -o ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_Bonf_signif_bins.bed \
#       # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
#       # echo -e "${pheno}.${CNV}.${filt}.perm_p" > ${list}
#       # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
#       # -o ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed \
#       # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
#       # rm ${list}
#       # #Make list of loci (±40kb merge distance & min 20kb size)
#       # ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
#       # bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct -d 40000 \
#       # -i ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed | \
#       # awk '{ if ($3-$2>20000) print $0 }' > \
#       # ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_nom_signif_bins.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_Bonf_signif_bins.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed
#     done
#   done
# done
# #Print table
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   for CNV in ANY_CNV CNV DEL DUP; do
#     for filt in ANY_FILTER all coding noncoding; do
#       if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz ]; then
#         echo -e "${pheno}\n${CNV}\n${filt}"
#         #Count total number of bins tested
#         zcat ${WRKDIR}/analysis/BIN_CNV_burdens/DD_vs_CTRL/DD_vs_CTRL_CNV_all.TBRden_results.bed.gz | \
#         sed '1d' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
#         #Count number of nominally significant bins from burden test
#         zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_nom_signif_bins.bed.gz | \
#         fgrep -v "#" | wc -l
#         #Count number of Bonferroni significant bins from burden test
#         zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.obs_Bonf_signif_bins.bed.gz | \
#         fgrep -v "#" | wc -l
#         #Count number of bins that passed permutation test
#         zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
#         fgrep -v "#" | wc -l
#         #Count number of merged loci that passed permutation test (as compared against master master master merged list of loci)
#         bedtools intersect -u \
#         -a <( zcat ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | cut -f1-3 ) \
#         -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
#       fi
#     done | paste - - - - - - - -
#   done
# done
# #Count overlap with syndromic CNV intervals
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   for CNV in ANY_CNV CNV DEL DUP; do
#     for filt in ANY_FILTER all coding noncoding; do
#       #Count raw overlaps with syndromic loci
#       bedtools intersect -u \
#       -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#       -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
#       bedtools intersect -wa -u -b - \
#       -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' ) | wc -l
#     done
#   done
# done | awk '{ print "=\""$1"/260\"" }'
# #Print overlaps with syndromic CNV intervals for non-coding CNVs
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   echo -e "\n${pheno}\n"
#   for CNV in ANY_CNV CNV DEL DUP; do
#       bedtools intersect -u \
#       -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#       -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
#       bedtools intersect -wa -u -b - \
#       -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' )
#   done | sort -Vk1,1 -k2,2n -k3,3n | uniq
# done

# #####Print top-5 most significant noncoding loci per comparison
# # for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
# for pheno in DD SCZ CNCR; do
#   echo -e "\n${pheno}\n"
#   while read chr start end; do
#     col=$( zcat ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz | \
#     head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${pheno}.ANY_CNV.noncoding.obs_p" ) | cut -f1 )
#     bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
#     sort -nk${col},${col} | head -n1 | cut -f1-4,${col}
#   done < <( bedtools intersect -u -wa \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
#     cut -f1-3 | fgrep -v "#" ) | sort -nrk5,5
# done

# #Prep for comparison of final loci within and between sets
# for pheno in DD SCZ DD_SCZ CNCR; do
#   #Make list of all loci
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
#     done
#   done | sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v pheno=${pheno} '{ print $0, pheno }' > ${TMPDIR}/${pheno}_master.bed
#   #Make lists of coding/noncoding loci for dels + dups
#   for filt in all coding noncoding; do
#     for CNV in DEL DUP; do
#       zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
#     done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/${pheno}_${filt}.bed
#   done
# done
# #Cat all master lists
# for pheno in DD SCZ DD_SCZ CNCR; do
#   cat ${TMPDIR}/${pheno}_master.bed
# done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/all_master.bed
# #Run comparison of final loci between del / dup within a given set
# for pheno in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       bedtools intersect -r -f 0.5 -c \
#       -a ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz \
#       -b ${TMPDIR}/${pheno}_${filt}.bed | awk '{ if ($NF==1) print $0 }' | wc -l
#     done
#   done
# done
# #Run comparison of loci specific to phenotype
# for pheno in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       awk -v pheno=${pheno} '{ if ($NF!=pheno) print $0 }' ${TMPDIR}/all_master.bed | \
#       bedtools intersect -v -r -f 0.5 -b - \
#       -a ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz | wc -l
#     done
#   done
# done

# #####Get list of most significant windows per locus per comparison (drawn vs master loci)
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   echo ${pheno}
#   for CNV in ANY_CNV CNV DEL DUP; do
#     echo ${CNV}
#     for filt in ANY_FILTER all coding noncoding; do
#       #Parallelize
#       bsub -q short -sla miket_sc -u nobody -J ${pheno}.${CNV}.${filt}.findPeakWindows \
#       "${WRKDIR}/bin/rCNVmap/bin/find_peak_window.sh ${pheno} ${CNV} ${filt}"
#       # while read chr start end; do
#       #   col=$( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
#       #   head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${pheno}.${CNV}.${filt}.obs_p" ) | cut -f1 )
#       #   bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
#       #   -a ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
#       #   sort -nk${col},${col} | head -n1 | cut -f1-4
#       # done < <( bedtools intersect -u -wa \
#       #   -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#       #   -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
#       #   cut -f1-3 | fgrep -v "#" ) > \
#       # ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
#     done
#   done
# done
# #Copy to directory for plotting
# if [ -e ${WRKDIR}/data/plot_data/signif_peak_windows ]; then
#   rm -rf ${WRKDIR}/data/plot_data/signif_peak_windows
# fi
# mkdir ${WRKDIR}/data/plot_data/signif_peak_windows
# for pheno in DD SCZ DD_SCZ CNCR ANY_DISEASE; do
#   for CNV in CNV DEL DUP ANY_CNV; do
#     for filt in all coding noncoding; do
#       cp ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed.gz \
#       ${WRKDIR}/data/plot_data/signif_peak_windows/
#     done
#   done
# done

# #####Get list of most representative windows per locus per comparison (drawn vs disease-specific loci)
# for pheno in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
#   echo ${pheno}
#   for CNV in ANY_CNV CNV DEL DUP; do
#     echo ${CNV}
#     for filt in ANY_FILTER all coding noncoding; do
#       #Parallelize
#       bsub -q short -sla miket_sc -u nobody -J ${pheno}.${CNV}.${filt}.findRepWindows \
#       "${WRKDIR}/bin/rCNVmap/bin/find_representative_windows.sh ${pheno} ${CNV} ${filt}"
#       # while read chr start end; do
#       #   size=$((${end}-${start}))
#       #   mod=$( expr ${size} % 100000 )
#       #   case ${mod} in
#       #     0)
#       #       paste <( seq ${start} 100000 $((${end}-100000)) ) \
#       #       <( seq $((${start}+100000)) 100000 ${end} ) | \
#       #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
#       #       bedtools intersect -wa -r -f 1 -b - \
#       #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
#       #       cut -f1-4
#       #       ;;
#       #     25000)
#       #       paste <( seq ${start} 100000 $((${end}-125000)) ) \
#       #       <( seq $((${start}+100000)) 100000 $((${end}-25000)) ) | \
#       #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
#       #       bedtools intersect -wa -r -f 1 -b - \
#       #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
#       #       cut -f1-4
#       #       ;;
#       #     50000)
#       #       paste <( seq $((${start}-25000)) 100000 $((${end}-75000)) ) \
#       #       <( seq $((${start}+75000)) 100000 $((${end}+25000)) ) | \
#       #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
#       #       bedtools intersect -wa -r -f 1 -b - \
#       #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
#       #       cut -f1-4
#       #       ;;
#       #     75000)
#       #       paste <( seq ${start} 100000 $((${end}-75000)) ) \
#       #       <( seq $((${start}+100000)) 100000 $((${end}+25000)) ) | \
#       #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
#       #       bedtools intersect -wa -r -f 1 -b - \
#       #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
#       #       cut -f1-4
#       #       ;;
#       #   esac
#       # done < <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
#       #   fgrep -v "#" | cut -f1-3 ) | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
#       # ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
#       # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
#     done
#   done
# done

# #####Get counts of significant loci per comparison for Figure 1d barplot
# for pheno in DD SCZ DD_SCZ CNCR; do
#   echo ${pheno}
#   for filt in ANY_FILTER coding noncoding; do
#     #Any CNV
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | cut -f1-3 )| wc -l
#     #DEL-specific
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
#     bedtools intersect -v -a - \
#     -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
#     #DUP-specific
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
#     bedtools intersect -v -a - \
#     -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
#   done | paste - - -
# done | paste - - - - > \
# ${WRKDIR}/data/plot_data/signif_loci_breakdown.counts.txt

# #####Get overlap for 3-way venn between CNCR, DD, SCZ
# for filt in ANY_FILTER coding noncoding; do
#   echo -e "\n${filt}"
#   for pheno in DD SCZ CNCR; do
#     echo "${pheno}"
#     notA=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${pheno} | head -n1 )
#     notB=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${pheno} | tail -n1 )
#     #Total
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
#     #Unique to pheno
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
#     #Shared with A but not B
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -u -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
#     #Shared with B but not A
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -u -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
#     #Shared with A and B
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -u -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
#     bedtools intersect -u -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
#   done | paste - - - - - -
# done

# #####Data for contrasts between deletion and duplication hotspots (Fig 1e)
# if [ -e ${WRKDIR}/data/plot_data/del_vs_dup ]; then
#   rm -rf ${WRKDIR}/data/plot_data/del_vs_dup
# fi
# mkdir ${WRKDIR}/data/plot_data/del_vs_dup
# #Get list of loci that are DEL only
# bedtools intersect -u \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
# bedtools intersect -v -a - \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz > \
# ${TMPDIR}/DEL_only.loci.bed
# bedtools intersect -u -b ${TMPDIR}/DEL_only.loci.bed \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
# ${TMPDIR}/DEL_only.peak_windows.bed
# #Get list of loci that are DUP only
# bedtools intersect -u \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
# bedtools intersect -v -a - \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
# ${TMPDIR}/DUP_only.loci.bed
# bedtools intersect -u -b ${TMPDIR}/DUP_only.loci.bed \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
# ${TMPDIR}/DUP_only.peak_windows.bed
# #Get list of loci that are DEL and DUP
# bedtools intersect -u \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
# bedtools intersect -u -a - \
# -b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
# ${TMPDIR}/DEL_and_DUP.loci.bed
# bedtools intersect -u -b ${TMPDIR}/DEL_and_DUP.loci.bed \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
# ${TMPDIR}/DEL_and_DUP.peak_windows.bed
# #Calculate average sizes
# for CNV in DEL_only DUP_only DEL_and_DUP; do
#   awk '{ print $3-$2 }' ${TMPDIR}/${CNV}.loci.bed > \
#   ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.size.txt
# done
# #Fraction of DEL-only and DUP-only hotspots per disease
# for pheno in CNCR DD_SCZ; do
#   for dummy in 1; do
#     echo ${pheno}
#     #ALL
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
#     #DEL only
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
#     #DUP only
#     bedtools intersect -u \
#     -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
#     bedtools intersect -v -a - \
#     -b ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
#   done | paste - - - -
# done > ${WRKDIR}/data/plot_data/del_vs_dup/del_vs_dup.count_by_class.txt
# #Genes per 100kb per hotspot
# for CNV in DEL_only DUP_only DEL_and_DUP; do
#   while read chr start end skip; do
#     size=$((${end}-${start}))
#     fgrep -wf ${SFARI_ANNO}/genelists/Gencode_proteinCoding.genes.list \
#     ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
#     bedtools intersect -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | \
#     awk '$4 !~ /\-AS1/ { print $4 }' | sort | uniq | wc -l | \
#     awk -v size=${size} '{ print (100000*$1)/size }'
#   done < ${TMPDIR}/${CNV}.loci.bed > \
#   ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.pc_genes_per_100kb.txt
# done
# #SD coverage per 100kb per hotspot
# for CNV in DEL_only DUP_only DEL_and_DUP; do
#   bedtools coverage -b <( cut -f1-3 ${TMPDIR}/${CNV}.loci.bed ) \
#   -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed ) | \
#   cut -f7 > ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.SD_coverage.txt
# done
# #Get max OR per hotspot per CNV class
# for CNV in DEL DUP; do
#   while read chr start end skip; do
#     for pheno in CNCR DD SCZ; do
#       for filt in all coding noncoding; do
#         bedtools intersect -wa \
#         -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | sed '1d' ) \
#         -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
#       done
#     done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
#   done < ${TMPDIR}/${CNV}_only.loci.bed > \
#   ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.peakORs.txt
# done
# while read chr start end skip; do
#   for pheno in CNCR DD SCZ; do
#     for filt in all coding noncoding; do
#       bedtools intersect -wa \
#       -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${pheno}_vs_CTRL/${pheno}_vs_CTRL_CNV_${filt}.TBRden_results.bed.gz | sed '1d' ) \
#       -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
#     done
#   done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
# done < ${TMPDIR}/DEL_and_DUP.loci.bed > \
# ${WRKDIR}/data/plot_data/del_vs_dup/DEL_and_DUP.peakORs.txt

# #####Count number of autosomal significant loci that don't overlap with consensus syndromic list
# for pheno in DD SCZ DD_SCZ; do
#   zcat ${WRKDIR}/analysis/Final_Loci/significant/${pheno}/${pheno}.ANY_CNV.ANY_FILTER.perm_signif_bins.bed.gz | \
#   fgrep -v "#" | cut -f1-3
# done | sort -Vk1,1 -k2,2n -k3,3n | uniq | bedtools intersect -u -wa -b - \
# -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | \
# cut -f1-3 | bedtools intersect -u -a - \
# -b ${WRKDIR}/data/unfiltered_annotations/CollinsFocal_PathogenicCNVs_allSources.boundaries.bed.gz | wc -l
