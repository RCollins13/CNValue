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
        bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 250000 -s 25000 -d 1000000 -r 0 \
        -o ${WRKDIR}/analysis/BIN_CNV_pileups/${pheno}/${pheno}.${CNV}.${VF}.${filt}.BIN_CNV_pileup.bed \
        -x ${WRKDIR}/data/master_annotations/other/hotspotAnalysis.excluded_loci.bed  \
        ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
        /data/talkowski/rlc47/src/GRCh37.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | cut -f1 )


#Run rCNV TBRden pileups (5kb bins, 5kb step, 1Mb flank dist; run at 0kb, 50kb, and 100kb smoothing)
while read group eti tier descrip include exclude color n; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${group} ]; then
    rm -r ${WRKDIR}/analysis/BIN_CNV_pileups/${group}
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${group}
  for smooth in 0 5 10; do
    if ! [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed ]; then
      mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed
    fi
    for CNV in DEL DUP CNV; do
      for filter in all coding dispensable noncoding intergenic; do
        #Parallelize intersections (LSF)
        bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 5000 -s 5000 -d 1000000 -r ${smooth} \
        -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.${filter}.bed \
        -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed \
        ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.GRCh37.${filter}.bed.gz \
        /data/talkowski/rlc47/src/GRCh37.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )
#Run urCNV TBRden pileups (50kb bins, 50kb step, 1Mb flank dist; run at 0kb, 100kb, and 250kb smoothing)
while read group eti tier descrip include exclude color n; do
  for smooth in 0 10 25; do
    if ! [ -e ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}kb_smoothed ]; then
      mkdir ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}kb_smoothed
    fi
    for CNV in DEL DUP CNV; do
      for filter in all coding dispensable noncoding intergenic; do
        #Parallelize intersections (LSF)
        bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_binned_pileup \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_binned_pileup.sh -z -w 50000 -s 50000 -d 1000000 -r ${smooth} \
        -o ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}kb_smoothed/${group}.${CNV}.urCNVs.TBRden_binned_pileup.${filter}.bed \
        -x ${WRKDIR}/lists/rCNVmap_excluded_loci.bins.bed \
        ${WRKDIR}/data/CNV/CNV_MASTER/${group}/${group}.${CNV}.urCNVs.GRCh37.${filter}.bed.gz \
        /data/talkowski/rlc47/src/GRCh37.genome"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )
#Rename output directories
while read group eti tier descrip include exclude color n; do
  for smooth in 5 10; do
    mv ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}0kb_smoothed
  done
  for smooth in 10 25; do
    mv ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}kb_smoothed \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}0kb_smoothed
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list )

#####Run TBRden analysis
#rCNVs
while read group eti tier descrip include exclude color n; do
  if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
    rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  fi
  mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  for smooth in 0 50 100; do
    if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed ]; then
      rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed
    fi
    mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed
    for CNV in DEL DUP CNV; do
      for filter in all coding dispensable noncoding intergenic; do
        #Parallelize analyses (LSF)
        bsub -q normal -sla miket_sc -J ${group}_${CNV}_TBRden_analysis -u nobody \
        -o ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.out \
        -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.err \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
        ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/${smooth}kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/ \
        ${group}_vs_CTRL_${CNV}_${filter} 0.0000003818455 ${color}"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" )
#urCNVs
while read group eti tier descrip include exclude color n; do
  for smooth in 0 100 250; do
    if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed ]; then
      rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed
    fi
    mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed
    for CNV in DEL DUP CNV; do
      for filter in all coding dispensable noncoding intergenic; do
        #Parallelize analyses (LSF)
        bsub -q short -sla miket_sc -J ${group}_${CNV}_TBRden_analysis -u nobody \
        -o ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.out \
        -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.err \
        "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
        ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/urCNV_${smooth}kb_smoothed/CTRL.${CNV}.urCNVs.TBRden_binned_pileup.${filter}.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/urCNV_${smooth}kb_smoothed/${group}.${CNV}.urCNVs.TBRden_binned_pileup.${filter}.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_${smooth}kb_smoothed/ \
        ${group}_vs_CTRL_${CNV}_${filter} 0.0000003818455 ${color}"
      done
    done
  done
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" )

#####Get significant urCNV loci from 100kb smoothed pileups, allowing for up to 50kb between bins that are merged
for CNV in CNV DEL DUP; do
  echo -e "\n\n${CNV}\n\n"
  while read group eti tier descrip include exclude color n; do
    echo -e "${group}\t${n}\t${descrip}"
    for filter in all coding dispensable noncoding intergenic; do
      if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_100kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.TBRden_results.bed.gz ]; then
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/urCNV_100kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.TBRden_results.bed.gz | \
        awk '{ if ($NF<=0.0000009545993) print $0 }' | bedtools merge -d 50000 -i - | wc -l
      else
        echo "NA"
      fi
    done
  done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
    fgrep -v "CTRL" ) | paste - - - - - -
done

# #####Run TBRden analysis (COE METHOD)
# #rCNVs
# while read group eti tier descrip include exclude color n; do
#   if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
#     rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
#   fi
#   mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
#   for smooth in 0 50 100; do
#     if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed ]; then
#       rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed
#     fi
#     mkdir ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed
#     for CNV in DEL DUP CNV; do
#       for filter in all coding dispensable noncoding intergenic; do
#         #Parallelize analyses (LSF)
#         bsub -q normal -sla miket_sc -J ${group}_${CNV}_TBRden_analysis -u nobody \
#         -o ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.out \
#         -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/${group}_vs_CTRL_${CNV}_${filter}.err \
#         "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
#         ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL/${smooth}kb_smoothed/CTRL.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
#         ${WRKDIR}/analysis/BIN_CNV_pileups/${group}/${smooth}kb_smoothed/${group}.${CNV}.TBRden_binned_pileup.${filter}.bed.gz \
#         ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${smooth}kb_smoothed/ \
#         ${group}_vs_CTRL_${CNV}_${filter} 0.0000003818455 ${color}"
#       done
#     done
#   done
# done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | fgrep -v "CTRL" )





















#####Run 1k CNV shift direct permutation tests for all comparisons
#Note: changed from old 100k matched Fisher permutation
#Note: initial p-value cutoff used: 0.05/52378 = 9.545993e-07 = 0.0000009545993
#This corresponds to the number of non-overlapping 50kb bins we tested (after blacklisting N-mask, etc)
#50kb chosen since minimum CNV size = 50kb, so max # independent tests = size of genome / 50kb
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
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_all \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
    #Coding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -c -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_coding \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
    #Noncoding CNVs only
    zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.TBRden_results.bed.gz | \
    awk -v OFS="\t" '{ if ($NF<=(0.05/130942.8) && $NF!="NA") print $0 }' | sort -Vk1,1 -k2,2n -k3,3n > \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed
    bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -z -n -N 1000 -t upper \
    -p ${group}_vs_CTRL_${CNV}_noncoding \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed \
    ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/"
  done
done
# #OLD CODE USED FOR FISHER PERMUTATION
#Split into partitions of 1k permutations each (x100 per comparison)
# for i in $( seq -w 001 100 ); do
#   echo ${i}
#   #Coding + noncoding CNVs
#   if ! [ -e ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i} ]; then
#     mkdir ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}
#   fi
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.all.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 1000000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_all.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_all \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_all.Bonferroni.bed"
#   #Coding CNVs only
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.coding.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_coding.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_coding \
#   -c -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_coding.Bonferroni.bed"
#   #Noncoding CNVs only
#   bsub -q short -sla miket_sc -J ${group}_vs_CTRL.${CNV}.noncoding.1k_permute.${i} -u nobody \
#   "${WRKDIR}/bin/rCNVmap/bin/CNV_shift_test.sh -z -d 5 -b 100000 -N 100 \
#   -o ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_noncoding.permuted.${i}.bed \
#   -p ${group}_vs_CTRL_${CNV}_noncoding \
#   -n -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
#   ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
#   ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_noncoding.Bonferroni.bed"
# done


#####Collect results from permutation tests
#Collect results across all permutations per group
# for group in DD SCZ DD_SCZ CNCR; do
#   for CNV in CNV DEL DUP; do
#     for filt in all coding noncoding; do
#       for i in $( seq -w 001 100 ); do
#         zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/${i}/${group}_vs_CTRL_${CNV}_${filt}.permuted.${i}.bed.gz | \
#         fgrep -v "#" | cut -f7 | paste -s
#       done > ${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt
#       Rscript -e "options(scipen=100);\
#       d <- apply(read.table(\"${TMPDIR}/${group}_${CNV}_${filt}_perm.in.txt\",header=F),2,sum);\
#       write.table(data.frame(\"perms_less_sig\"=10000-d,\"perms_as_or_more_sig\"=d,\"perm_p\"=(d+1)/10001),\
#         \"${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt\",row.names=F,col.names=T,sep=\"\t\",quote=F)"
#       paste <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/perm_split/001/${group}_vs_CTRL_${CNV}_${filt}.permuted.001.bed.gz | cut -f1-5 ) \
#       ${TMPDIR}/${group}_${CNV}_${filt}_perm.out.txt > \
#       ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.permuted.merged.bed
#     done
#   done
# done
#Compose table of all bins for each comparison
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      cat <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -wb -f 1 -r -a - \
      -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
      cut -f1-17,26 ) \
      <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
      sed '1d' | bedtools intersect -v -f 1 -r -a - \
      -b <( zcat ${WRKDIR}/analysis/BIN_CNV_permutation/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_direct_test_results.bed.gz | sed '1d' ) | \
      awk -v OFS="\t" '{ print $0, "NA" }' ) | sort -Vk1,1 -k2,2n -k3,3n | \
      cat <( paste <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | \
        head -n1 | awk '{ print "#"$0 }' ) <( echo -e "perm_sim_p" ) ) - > \
      ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
      gzip -f ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed
    done
  done
done
#Make master table of all bins with summary stats for each bin & each comparison
for dummy in 1; do
  for group in CTRL DD SCZ CNCR; do
    for CNV in DEL DUP; do
      zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
      fgrep -v "#" | cut -f5 > ${TMPDIR}/${group}.${CNV}.all.CNVcounts.tmp
      for filt in coding noncoding; do
        zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.TBRden_binned_pileup.${filt}.bed.gz | \
        fgrep -v "#" | cut -f5 > ${TMPDIR}/${group}.${CNV}.${filt}.CNVcounts.tmp
        echo ${TMPDIR}/${group}.${CNV}.${filt}.CNVcounts.tmp
      done
    done
  done
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all coding noncoding; do
        zcat ${WRKDIR}/analysis/Final_Loci/${group}_vs_CTRL_${CNV}_${filt}.results.all_bins.bed.gz | \
        fgrep -v "#" | cut -f17-18 > ${TMPDIR}/${group}.${CNV}.${filt}.pvals.tmp
        echo ${TMPDIR}/${group}.${CNV}.${filt}.pvals.tmp
      done
    done
  done
done > ${TMPDIR}/list_of_paths.txt
for dummy in 1; do
  zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
  cut -f1-4 | head -n1
  for group in CTRL DD SCZ CNCR; do
    for CNV in DEL DUP; do
      for filt in coding noncoding; do
        echo "${group}.${CNV}.${filt}.count"
      done
    done
  done
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in CNV DEL DUP; do
      for filt in all coding noncoding; do
        echo -e "${group}.${CNV}.${filt}.obs_p\n${group}.${CNV}.${filt}.perm_p"
      done
    done
  done
done | paste -s > ${TMPDIR}/header.txt
paste <( zcat ${WRKDIR}/analysis/Final_Loci/DD_vs_CTRL_CNV_all.results.all_bins.bed.gz | \
  fgrep -v "#" | cut -f1-4 ) $( cat ${TMPDIR}/list_of_paths.txt ) | \
cat ${TMPDIR}/header.txt - > ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
${WRKDIR}/bin/rCNVmap/bin/summarize_master_burden_file.R \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2
mv ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed2 \
${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed
gzip -f ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed

#####Cut master table of bins that are significant per disease 
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${group} ]; then
    rm -rf ${WRKDIR}/analysis/Final_Loci/significant/${group}
  fi
  mkdir ${WRKDIR}/analysis/Final_Loci/significant/${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      echo ${filt}
      #Parallelize:
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_filterMasterBurden \
      "${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file_parallelized.sh \
      ${group} ${CNV} ${filt}"
      # list=`mktemp`
      # echo -e "${group}.${CNV}.${filt}.obs_p" > ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R -t 0.05 \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # echo -e "${group}.${CNV}.${filt}.perm_p" > ${list}
      # ${WRKDIR}/bin/rCNVmap/bin/filter_master_burden_file.R \
      # -o ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed \
      # ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz ${list}
      # rm ${list}
      # #Make list of loci (±40kb merge distance & min 20kb size)
      # ncol=$( head -n1 ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | awk '{ print NF }' )
      # bedtools merge -header -c $( seq 4 ${ncol} | paste -s -d, ) -o distinct -d 40000 \
      # -i ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed | \
      # awk '{ if ($3-$2>20000) print $0 }' > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed
    done
  done
done
#Print table
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    for filt in ANY_FILTER all coding noncoding; do
      if [ -e ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz ]; then
        echo -e "${group}\n${CNV}\n${filt}"
        #Count total number of bins tested
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/DD_vs_CTRL/DD_vs_CTRL_CNV_all.TBRden_results.bed.gz | \
        sed '1d' | awk '{ if ($1!="X" && $1!="Y") print $0 }' | wc -l
        #Count number of nominally significant bins from burden test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_nom_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of Bonferroni significant bins from burden test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.obs_Bonf_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of bins that passed permutation test
        zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
        fgrep -v "#" | wc -l
        #Count number of merged loci that passed permutation test (as compared against master master master merged list of loci)
        bedtools intersect -u \
        -a <( zcat ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | cut -f1-3 ) \
        -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
      fi
    done | paste - - - - - - - -
  done
done
#Count overlap with syndromic CNV intervals
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  for CNV in ANY_CNV CNV DEL DUP; do
    for filt in ANY_FILTER all coding noncoding; do
      #Count raw overlaps with syndromic loci
      bedtools intersect -u \
      -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
      bedtools intersect -wa -u -b - \
      -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' ) | wc -l
    done
  done
done | awk '{ print "=\""$1"/260\"" }'
#Print overlaps with syndromic CNV intervals for non-coding CNVs
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo -e "\n${group}\n"
  for CNV in ANY_CNV CNV DEL DUP; do
      bedtools intersect -u \
      -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      bedtools intersect -wa -u -b - \
      -a <( cut -f1,4- ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19.bed | sed '1d' )
  done | sort -Vk1,1 -k2,2n -k3,3n | uniq
done

#####Print top-5 most significant noncoding loci per comparison
# for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
for group in DD SCZ CNCR; do
  echo -e "\n${group}\n"
  while read chr start end; do
    col=$( zcat ${WRKDIR}/analysis/Final_Loci/MASTER.p_values.all_bins.bed.gz | \
    head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${group}.ANY_CNV.noncoding.obs_p" ) | cut -f1 )
    bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
    -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
    sort -nk${col},${col} | head -n1 | cut -f1-4,${col}
  done < <( bedtools intersect -u -wa \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.noncoding.perm_signif_bins.bed.gz | \
    cut -f1-3 | fgrep -v "#" ) | sort -nrk5,5
done

#Prep for comparison of final loci within and between sets
for group in DD SCZ DD_SCZ CNCR; do
  #Make list of all loci
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
    done
  done | sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" -v group=${group} '{ print $0, group }' > ${TMPDIR}/${group}_master.bed
  #Make lists of coding/noncoding loci for dels + dups
  for filt in all coding noncoding; do
    for CNV in DEL DUP; do
      zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | fgrep -v "#"
    done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/${group}_${filt}.bed
  done
done
#Cat all master lists
for group in DD SCZ DD_SCZ CNCR; do
  cat ${TMPDIR}/${group}_master.bed
done | sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/all_master.bed
#Run comparison of final loci between del / dup within a given set
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      bedtools intersect -r -f 0.5 -c \
      -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
      -b ${TMPDIR}/${group}_${filt}.bed | awk '{ if ($NF==1) print $0 }' | wc -l
    done
  done
done
#Run comparison of loci specific to phenotype
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    for filt in all coding noncoding; do
      awk -v group=${group} '{ if ($NF!=group) print $0 }' ${TMPDIR}/all_master.bed | \
      bedtools intersect -v -r -f 0.5 -b - \
      -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | wc -l
    done
  done
done

#####Get list of most significant windows per locus per comparison (drawn vs master loci)
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      #Parallelize
      bsub -q short -sla miket_sc -u nobody -J ${group}.${CNV}.${filt}.findPeakWindows \
      "${WRKDIR}/bin/rCNVmap/bin/find_peak_window.sh ${group} ${CNV} ${filt}"
      # while read chr start end; do
      #   col=$( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
      #   head -n1 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w $( echo -e "${group}.${CNV}.${filt}.obs_p" ) | cut -f1 )
      #   bedtools intersect -wa -u -b <( echo -e "${chr}\t${start}\t${end}" ) \
      #   -a ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      #   sort -nk${col},${col} | head -n1 | cut -f1-4
      # done < <( bedtools intersect -u -wa \
      #   -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
      #   -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_bins.bed.gz | \
      #   cut -f1-3 | fgrep -v "#" ) > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed
    done
  done
done
#Copy to directory for plotting
if [ -e ${WRKDIR}/data/plot_data/signif_peak_windows ]; then
  rm -rf ${WRKDIR}/data/plot_data/signif_peak_windows
fi
mkdir ${WRKDIR}/data/plot_data/signif_peak_windows
for group in DD SCZ DD_SCZ CNCR ANY_DISEASE; do
  for CNV in CNV DEL DUP ANY_CNV; do
    for filt in all coding noncoding; do
      cp ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.peak_windows.bed.gz \
      ${WRKDIR}/data/plot_data/signif_peak_windows/
    done
  done
done

#####Get list of most representative windows per locus per comparison (drawn vs disease-specific loci)
for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for CNV in ANY_CNV CNV DEL DUP; do
    echo ${CNV}
    for filt in ANY_FILTER all coding noncoding; do
      #Parallelize
      bsub -q short -sla miket_sc -u nobody -J ${group}.${CNV}.${filt}.findRepWindows \
      "${WRKDIR}/bin/rCNVmap/bin/find_representative_windows.sh ${group} ${CNV} ${filt}"
      # while read chr start end; do
      #   size=$((${end}-${start}))
      #   mod=$( expr ${size} % 100000 )
      #   case ${mod} in
      #     0)
      #       paste <( seq ${start} 100000 $((${end}-100000)) ) \
      #       <( seq $((${start}+100000)) 100000 ${end} ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     25000)
      #       paste <( seq ${start} 100000 $((${end}-125000)) ) \
      #       <( seq $((${start}+100000)) 100000 $((${end}-25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     50000)
      #       paste <( seq $((${start}-25000)) 100000 $((${end}-75000)) ) \
      #       <( seq $((${start}+75000)) 100000 $((${end}+25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #     75000)
      #       paste <( seq ${start} 100000 $((${end}-75000)) ) \
      #       <( seq $((${start}+100000)) 100000 $((${end}+25000)) ) | \
      #       awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }' | \
      #       bedtools intersect -wa -r -f 1 -b - \
      #       -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
      #       cut -f1-4
      #       ;;
      #   esac
      # done < <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz | \
      #   fgrep -v "#" | cut -f1-3 ) | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
      # ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
      # gzip -f ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.representative_windows.bed
    done
  done
done

#####Get counts of significant loci per comparison for Figure 1d barplot
for group in DD SCZ DD_SCZ CNCR; do
  echo ${group}
  for filt in ANY_FILTER coding noncoding; do
    #Any CNV
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | cut -f1-3 )| wc -l
    #DEL-specific
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
    bedtools intersect -v -a - \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
    #DUP-specific
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | \
    bedtools intersect -v -a - \
    -b <( zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.${filt}.perm_signif_bins.bed.gz | cut -f1-3 ) | wc -l
  done | paste - - -
done | paste - - - - > \
${WRKDIR}/data/plot_data/signif_loci_breakdown.counts.txt

#####Get overlap for 3-way venn between CNCR, DD, SCZ
for filt in ANY_FILTER coding noncoding; do
  echo -e "\n${filt}"
  for group in DD SCZ CNCR; do
    echo "${group}"
    notA=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${group} | head -n1 )
    notB=$( echo -e "DD\nSCZ\nCNCR" | fgrep -v ${group} | tail -n1 )
    #Total
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Unique to group
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with A but not B
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with B but not A
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
    #Shared with A and B
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notA}/${notA}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | \
    bedtools intersect -u -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${notB}/${notB}.ANY_CNV.${filt}.perm_signif_bins.bed.gz | wc -l
  done | paste - - - - - -
done

#####Data for contrasts between deletion and duplication hotspots (Fig 1e)
if [ -e ${WRKDIR}/data/plot_data/del_vs_dup ]; then
  rm -rf ${WRKDIR}/data/plot_data/del_vs_dup
fi
mkdir ${WRKDIR}/data/plot_data/del_vs_dup
#Get list of loci that are DEL only
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -v -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DEL_only.loci.bed
bedtools intersect -u -b ${TMPDIR}/DEL_only.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DEL_only.peak_windows.bed
#Get list of loci that are DUP only
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -v -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DUP_only.loci.bed
bedtools intersect -u -b ${TMPDIR}/DUP_only.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DUP_only.peak_windows.bed
#Get list of loci that are DEL and DUP
bedtools intersect -u \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
bedtools intersect -u -a - \
-b ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.DEL.ANY_FILTER.perm_signif_loci.bed.gz > \
${TMPDIR}/DEL_and_DUP.loci.bed
bedtools intersect -u -b ${TMPDIR}/DEL_and_DUP.loci.bed \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.peak_windows.bed.gz > \
${TMPDIR}/DEL_and_DUP.peak_windows.bed
#Calculate average sizes
for CNV in DEL_only DUP_only DEL_and_DUP; do
  awk '{ print $3-$2 }' ${TMPDIR}/${CNV}.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.size.txt
done
#Fraction of DEL-only and DUP-only hotspots per disease
for group in CNCR DD_SCZ; do
  for dummy in 1; do
    echo ${group}
    #ALL
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
    #DEL only
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
    #DUP only
    bedtools intersect -u \
    -a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DUP.ANY_FILTER.perm_signif_loci.bed.gz | \
    bedtools intersect -v -a - \
    -b ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.DEL.ANY_FILTER.perm_signif_loci.bed.gz | wc -l
  done | paste - - - -
done > ${WRKDIR}/data/plot_data/del_vs_dup/del_vs_dup.count_by_class.txt
#Genes per 100kb per hotspot
for CNV in DEL_only DUP_only DEL_and_DUP; do
  while read chr start end skip; do
    size=$((${end}-${start}))
    fgrep -wf ${SFARI_ANNO}/genelists/Gencode_proteinCoding.genes.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
    bedtools intersect -a - -b <( echo -e "${chr}\t${start}\t${end}" ) | \
    awk '$4 !~ /\-AS1/ { print $4 }' | sort | uniq | wc -l | \
    awk -v size=${size} '{ print (100000*$1)/size }'
  done < ${TMPDIR}/${CNV}.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.pc_genes_per_100kb.txt
done
#SD coverage per 100kb per hotspot
for CNV in DEL_only DUP_only DEL_and_DUP; do
  bedtools coverage -b <( cut -f1-3 ${TMPDIR}/${CNV}.loci.bed ) \
  -a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed ) | \
  cut -f7 > ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.SD_coverage.txt
done
#Get max OR per hotspot per CNV class
for CNV in DEL DUP; do
  while read chr start end skip; do
    for group in CNCR DD SCZ; do
      for filt in all coding noncoding; do
        bedtools intersect -wa \
        -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.TBRden_results.bed.gz | sed '1d' ) \
        -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
      done
    done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
  done < ${TMPDIR}/${CNV}_only.loci.bed > \
  ${WRKDIR}/data/plot_data/del_vs_dup/${CNV}.peakORs.txt
done
while read chr start end skip; do
  for group in CNCR DD SCZ; do
    for filt in all coding noncoding; do
      bedtools intersect -wa \
      -a <( zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_CNV_${filt}.TBRden_results.bed.gz | sed '1d' ) \
      -b <( echo -e "${chr}\t${start}\t${end}" ) | awk -v OFS="\n" '{ print $12, $14, $16 }'
    done
  done | sort -Vrk1,1 | fgrep -v Inf | fgrep -v NA | head -n1
done < ${TMPDIR}/DEL_and_DUP.loci.bed > \
${WRKDIR}/data/plot_data/del_vs_dup/DEL_and_DUP.peakORs.txt

#####Count number of autosomal significant loci that don't overlap with consensus syndromic list
for group in DD SCZ DD_SCZ; do
  zcat ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.ANY_CNV.ANY_FILTER.perm_signif_bins.bed.gz | \
  fgrep -v "#" | cut -f1-3
done | sort -Vk1,1 -k2,2n -k3,3n | uniq | bedtools intersect -u -wa -b - \
-a ${WRKDIR}/analysis/Final_Loci/significant/ANY_DISEASE/ANY_DISEASE.ANY_CNV.ANY_FILTER.perm_signif_loci.bed.gz | \
cut -f1-3 | bedtools intersect -u -a - \
-b ${WRKDIR}/data/unfiltered_annotations/CollinsFocal_PathogenicCNVs_allSources.boundaries.bed.gz | wc -l

#####Prepare unfiltered annotation files for enrichment tests
#SFARI gene sets
while read genes; do
  echo ${genes}
  fgrep -wf ${SFARI_ANNO}/genelists/${genes}.genes.list \
  ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | \
  awk '$4 !~ /\-AS/ { print $0 }' | grep -ve '^GL\|^X\|^Y\|^M' | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/${genes}.boundaries.bed
done < <( l ${SFARI_ANNO}/genelists/*genes.list | awk '{ print $9 }' | sed 's/\.genes\.list//g' | \
  xargs -I {} basename {} | fgrep -wv genes_merged_list | fgrep -wv Landrum2014_clinvar_XLinkedDiseaseAssociated )
for group in ASD EPI ID; do
  for level in HC MC LC; do
    echo -e "${group}_${level}"
    #Count
    fgrep -wf /data/talkowski/Samples/SFARI/EcoFinal/annotation_files/genes/${group}_${level}.list \
    ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed |  \
    awk '$4 !~ /\-AS/ { print $0 }' | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
    ${WRKDIR}/data/unfiltered_annotations/${genes}.boundaries.bed
  done
done
#GTEx top 1k most highly expressed genes by tissue
while read tissue; do
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed.gz | \
  sort -nrk5,5 | grep -ve '^GL\|^X\|^Y\|^M' | head -n1000 | cut -f4 | sort | uniq | \
  fgrep -wf - ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed |  awk '$4 !~ /\-AS/ { print $0 }' | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_highExpressors_top1k_${tissue}.boundaries.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
#UCNEs
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Dimitrieva2013_UCNE.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/UCNEs.boundaries.bed
#RepeatMasker by class
for class in LINE SINE LTR Satellite Simple_repeat Low_complexity; do
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/src/GRCh37.${class}.RMSK.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/RepMask_${class}.boundaries.bed
done
#HARs
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Doan2016_HAR.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/HAR.boundaries.bed
#SegDups
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/SegmentalDuplications_GRCh37.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/SegDups.boundaries.bed
#ClinGen Pathogenic CNV Regions
zcat ${SFARI_ANNO}/misc/ClinGen_PathogenicCNVs.bed.gz | \
grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/ClinGen_PathoCNV.boundaries.bed
#GWAS catalog variants
zcat ${SFARI_ANNO}/misc/GWAS_Catalog_variants.merged_simple.bed.gz | \
grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/GWAS_Catalog_variants.boundaries.bed
#ENCODE TF ChIP peaks
while read TF; do
  echo ${TF}
  for level in all 0_250 250_500 500_750 750_1000; do
    zcat ${SFARI_ANNO}/noncoding/ENCODE/split_TF_peaks/${TF}/${TF}.ENCODE_peaks.${level}.bed.gz | \
    grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
    ${WRKDIR}/data/unfiltered_annotations/ENCODE_${TF}_peaks_${level}.boundaries.bed
  done
done < <( zcat ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeRegTfbsClusteredV3.bed.gz | cut -f4 | sort | uniq | \
  cat <( echo -e "all_TF" ) - )
#GTEx eQTLs
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_eQTLs_${tissue}.boundaries.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/GTEx_eQTLs/${tissue}/${tissue}.signif_eQTLs.bed
done < <( l ${SFARI_ANNO}/misc/GTEx_eQTLs/*snpgenes | awk '{ print $9 }' | \
  sed -e 's/GTEx_eQTLs\//\t/g' -e 's/_Analysis\.snpgenes//g' | awk '{ print $2 }' ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > \
  ${WRKDIR}/data/unfiltered_annotations/GTEx_eQTLs_AllTissues.boundaries.bed
#COSMIC variants
zcat ${SFARI_ANNO}/misc/COSMIC_variants.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/COSMIC_variants.boundaries.bed
#Affy6 Probes
zcat ${SFARI_ANNO}/misc/Affy6_Probes.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/Affy6_Probes.boundaries.bed
#Illumina 1MDuo Probes
zcat ${SFARI_ANNO}/misc/Illumina_1MDuo_Probes.bed.gz | grep -ve '^GL\|^X\|^Y\|^M' | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/Illumina_1MDuo_Probes.boundaries.bed
#Schmitt TBRs
while read tissue; do
  echo -e "${tissue}"
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/TBRs_${tissue}.boundaries.bed
done < <( cat <( echo "MERGED" ) /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' /data/talkowski/rlc47/TAD_intolerance/data/nuc/TBR/${tissue}.TBR.bed
done < /data/talkowski/rlc47/TAD_intolerance/lists/Schmitt_tissues.list | sort -Vk1,1 -k2,2n -k3,3n | uniq > \
${WRKDIR}/data/unfiltered_annotations/TBRs_AllTissues.boundaries.bed
#Super enhancers
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/superEnhancers/cleaned/${tissue}_superEnhancer.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/SuperEnhancers_${tissue}.boundaries.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/superEnhancers/*bed | \
  awk '{ print $9 }' | sed 's/superEnhancers\//\t/g' | sed 's/\.bed/\t/g' | cut -f2 )
#Tissue-specific enhancers
while read tissue; do
  echo ${tissue}
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Enhancers_${tissue}.boundaries.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 )
while read tissue; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/TissueEnhancers/cleaned/${tissue}_enhancers.bed
done < <( ls -l ${SFARI_ANNO}/noncoding/TissueEnhancers/*txt | \
  awk '{ print $9 }' | sed 's/TissueEnhancers\//\t/g' | sed 's/_EP/\t/g' | cut -f2 ) | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Enhancers_all.boundaries.bed
#Brain enhancer RNAs (count)
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Yao2015_BrainEnhancerRNA.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Brain_eRNAs.boundaries.bed
#FMRP targets
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Ascano2012_FMRP.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/FMRP_targets.boundaries.bed
#FANTOM5 NPC, fetal brain, and adult brain enhancers
for tissue in NPC AdultBrain FetalBrain; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/FANTOM52015_${tissue}_TPM5.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/FANTOM5_${tissue}_enhancers.boundaries.bed
done
#Cotney CHD8 peaks
for tissue in Brain hNSC; do
  grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/noncoding/Cotney2015_CHD8_${tissue}_peaks.bed | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/CotneyCHD8_${tissue}_peaks.boundaries.bed
done
#Talkowski 36 rCNV loci (count)
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/rCNVs/talkowski_highQual_rCNVs.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Talkowski_rCNVs_HQ.boundaries.bed
#Claire's new rCNV list
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_CER.bed | \
cut -f1,4-5 | sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/Redin_PathogenicCNVs_allSources.boundaries.bed
#Final consensus rCNV list
fgrep -v "#" ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $4, "CNV" }' | \
sort -Vk1,1 -k2,2n -k3,3n > ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed
bedtools intersect -r -f 0.5 -wa -wb \
-a ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed \
-b ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed > \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.selfIntersect.bed
cut -f4 ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.bed > ${TMPDIR}/CNV_interval_IDs.list
source /apps/lab/miket/anaconda/4.0.5/envs/collins_py3/bin/activate collins_py3
/data/talkowski/rlc47/code/svcf/scripts/bedcluster -p PathogenicCNVs_allSources_nonredundant_hg19_RLC -m \
${TMPDIR}/CNV_interval_IDs.list \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.selfIntersect.bed \
${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed
while read CNVID; do
  awk -v OFS="\t" -v CNVID=${CNVID} '{ if ($7==CNVID) print $1, $2, $3, $5 }' \
  ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed | \
  bedtools merge -c 4 -o distinct -i - 
done < <( fgrep -v "#" ${TMPDIR}/PathogenicCNVs_allSources_hg19_RLC.merged.bed | \
  cut -f7 | sort -Vk1,1 | uniq ) | sort -Vk1,1 -k2,2n -k3,3n > \
${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_RLC.bed
grep -ve '^GL\|^X\|^Y\|^M' ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_nonredundant_hg19_RLC.bed | \
sort -Vk1,1 -k2,2n -k3,3n | uniq > ${WRKDIR}/data/unfiltered_annotations/CollinsFocal_PathogenicCNVs_allSources.boundaries.bed
#rCNV list by source
while read study; do
  echo ${study}
  awk -v study=${study} '{ if ($4==study) print $0 }' \
  ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
  sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/unfiltered_annotations/${study}_PathoCNVs.boundaries.bed
done < <( fgrep -v "#" ${SFARI_ANNO}/misc/PathogenicCNVs_allSources_hg19_RLC.bed | \
  cut -f4 | sort | uniq )
#Gzip all annotations
gzip -f ${WRKDIR}/data/unfiltered_annotations/*.boundaries.bed

#####Launch unfiltered hotspot functional enrichment tests (both binomial and t-test)
if ! [ -e ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered ]; then
  mkdir ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered
fi
while read anno; do
  if [ -e ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/ ]; then
    rm -rf ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/
  fi
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/
  for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
    for CNV in ANY_CNV CNV DEL DUP; do
      for filt in ANY_FILTER all coding noncoding; do 
        #Binomial test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t binomial \
        -p ${group}_${CNV}_${filt}.${anno}.unfiltered_binomial \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/"
        #t test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t t \
        -p ${group}_${CNV}_${filt}.${anno}.unfiltered_t \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/"
      done
    done
  done > ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh
  chmod a+x ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh
  bsub -q short -sla miket_sc -J hotspot_unfiltered_enrichment_${anno} -u nobody \
  "sh ${WRKDIR}/bin/LSF/hotspot_enrichments_unfiltered/${anno}.submit.sh"
# done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
#         awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )
#optional: exclude ENCODE, GTEx, Enhancer, Super Enhancer, and TAD datasets
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve 'GTEx\|ENCODE\|Enhancer\|TBRs_' )

#####Launch raw CNV set enrichment tests
paste <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve '250_500\|500_750' ) \
<( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | grep -ve '250_500\|500_750' ) > \
${TMPDIR}/annotations_to_test.list
#Run first pass with adjustment to size
for group in DD SCZ DD_SCZ CNCR; do
  if [ -e ${WRKDIR}/analysis/CNV_set_enrichments/${group} ]; then
    rm -rf ${WRKDIR}/analysis/CNV_set_enrichments/${group}
  fi
  mkdir ${WRKDIR}/analysis/CNV_set_enrichments/${group}
  for CNV in CNV DEL DUP; do
    #All
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt"
    for filt in coding noncoding; do 
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt"
    done
  done
done
#Run second pass with adjustment to all genes, not to size
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    #All
    adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt.gz | \
      fgrep -w gencode.GRCh37.basic_gene_symbols | cut -f6 )
    echo ${adjust}
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment_geneNorm -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.geneNorm.setEnrichments.txt"
    for filt in coding noncoding; do 
      adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt.gz | \
        fgrep -w gencode.GRCh37.basic_gene_symbols | cut -f6 )
      echo ${adjust}
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment_geneNorm -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.geneNorm.setEnrichments.txt"
    done
  done
done
#Run third pass with adjustment based on median enrichment for all annotations tested
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in CNV DEL DUP; do
    #All
    adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.setEnrichments.txt.gz | \
      cut -f6 | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
    echo ${adjust}
    bsub -q short -sla miket_sc -J ${group}_${CNV}_all_setEnrichment_medEnrichNorm -u nobody \
    "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${TMPDIR}/annotations_to_test.list \
    ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_all.medEnrichNorm.setEnrichments.txt"
    for filt in coding noncoding; do 
      adjust=$( zcat ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.setEnrichments.txt.gz | \
        cut -f6 | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
      echo ${adjust}
      bsub -q short -sla miket_sc -J ${group}_${CNV}_${filt}_setEnrichment_medEnrichNorm -u nobody \
      "${WRKDIR}/bin/rCNVmap/bin/CNV_set_annoClass_bulk_test.sh -z -m ${adjust} \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.${filt}.bed.gz \
      ${TMPDIR}/annotations_to_test.list \
      ${WRKDIR}/analysis/CNV_set_enrichments/${group}/${group}_${CNV}_${filt}.medEnrichNorm.setEnrichments.txt"
    done
  done
done


#####Filter all annotations by excluding elements within 10kb of a protein-coding exon
while read anno; do
  echo ${anno}
  bedtools intersect -v -wa \
  -a ${WRKDIR}/data/unfiltered_annotations/${anno}.boundaries.bed.gz \
  -b ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.10kb_flanks.merged.bed > \
  ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed
  if [ -s ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed ]; then
    gzip -f ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed 
  else
    echo "Removing zero-sized annotation file for ${anno}"
    rm ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed
  fi
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | xargs -I {} basename {} | sed 's/\.boundaries\.bed\.gz//g' )

#####Launch filtered hotspot functional enrichment tests (both binomial and t-test) for noncoding sites only
if ! [ -e ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered ]; then
  mkdir ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered
fi
while read anno; do
  if [ -e ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/ ]; then
    rm -rf ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/
  fi
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/
  for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
    for CNV in ANY_CNV CNV DEL DUP; do
      for filt in noncoding; do 
        #Binomial test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t binomial \
        -p ${group}_${CNV}_${filt}.${anno}.filtered_binomial \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/"
        #t test
        echo -e "${WRKDIR}/bin/rCNVmap/bin/hotspot_annotation_test.sh -N 1000 -a greater -t t \
        -p ${group}_${CNV}_${filt}.${anno}.filtered_t \
        ${WRKDIR}/analysis/Final_Loci/significant/${group}/${group}.${CNV}.${filt}.perm_signif_loci.bed.gz \
        ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz \
        ${WRKDIR}/data/filtered_annotations/${anno}.boundaries.bed.gz \
        ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/"
      done
    done
  done > ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh
  chmod a+x ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh
  bsub -q short -sla miket_sc -J hotspot_filtered_enrichment_${anno} -u nobody \
  "sh ${WRKDIR}/bin/LSF/hotspot_enrichments_filtered/${anno}.submit.sh"
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' | \
        grep -ve '250_500\|500_750' )

#####Collect results from unfiltered enrichment tests
if ! [ -e ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS ]; then
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS
fi
#Print headers to output files
for dummy in 1; do
  for prefix in binomial_OR binomial_CI_min binomial_CI_max binomial_p \
  t_fold t_CI_min t_CI_max t_p; do
    for dummy in 2; do
      echo "annotation"
      for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
        for CNV in ANY_CNV CNV DEL DUP; do
          for filt in ANY_FILTER all coding noncoding; do
            echo -e "${group}_${CNV}_${filt}"
          done
        done
      done | paste -s
    done | paste - - > ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.${prefix}.txt
  done
done
while read anno; do
  echo ${anno}
  #Binomial odds ratio
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_OR.txt
  #Binomial odds ratio CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_CI_min.txt
  #Binomial odds ratio CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_CI_max.txt
  #Binomial p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.binomial_p.txt
  #t-test fold-enrichment
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_fold.txt
  #t-test CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_CI_min.txt
  #t-test CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_CI_max.txt
  #t-test p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/RESULTS/hotspot_unfiltered_annotations.t_p.txt
done < <( l ${WRKDIR}/data/unfiltered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )

#####Collect results from filtered enrichment tests
if ! [ -e ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS ]; then
  mkdir ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS
fi
#Print headers to output files
for dummy in 1; do
  for prefix in binomial_OR binomial_CI_min binomial_CI_max binomial_p \
  t_fold t_CI_min t_CI_max t_p; do
    for dummy in 2; do
      echo "annotation"
      for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
        for CNV in ANY_CNV CNV DEL DUP; do
          for filt in ANY_FILTER all coding noncoding; do
            echo -e "${group}_${CNV}_${filt}"
          done
        done
      done | paste -s
    done | paste - - > ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.${prefix}.txt
  done
done
while read anno; do
  echo ${anno}
  #Binomial odds ratio
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_OR.txt
  #Binomial odds ratio CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_CI_min.txt
  #Binomial odds ratio CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_CI_max.txt
  #Binomial p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_binomial.TBRden_binomial_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.binomial_p.txt
  #t-test fold-enrichment
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f4 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_fold.txt
  #t-test CI min
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f5 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_CI_min.txt
  #t-test CI max
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f6 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_CI_max.txt
  #t-test p value
  for dummy in 1; do
    echo ${anno}
    for group in ANY_DISEASE DD SCZ DD_SCZ CNCR; do
      for CNV in ANY_CNV CNV DEL DUP; do
        for filt in ANY_FILTER all coding noncoding; do
          outfile=${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/${anno}/${group}_${CNV}_${filt}.${anno}.filtered_t.TBRden_ttest_annotation_test.results.txt
          if [ -e ${outfile} ]; then
            cut -f7 ${outfile} | tail -n1
          else
            echo "NA"
          fi
        done | paste -s
      done | paste -s
    done | paste -s
  done | paste - - >> ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/RESULTS/hotspot_filtered_annotations.t_p.txt
done < <( l ${WRKDIR}/data/filtered_annotations/*boundaries.bed.gz | \
        awk '{ print $9 }' | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/\.boundaries\.bed\.gz//g' )

#####Gather data for published CNV locus enrichments (Fig1f)
#Get ORs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f4
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f4
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.ORs.txt
#Get lower CIs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.lowerCI.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f5
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f5
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.lowerCI.txt
#Get upper CIs
for dummy in 1; do
  echo "group"
  for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs \
  DDD_2016_PathoCNVs ACMG_2013_PathoCNVs \
  Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
    echo "${anno}"
  done
done | paste -s > ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.upperCI.txt
for group in ANY_DISEASE CNCR DD SCZ; do
  for dummy in 1; do
    echo ${group}
    for anno in CollinsFocal_PathogenicCNVs_allSources ClinGen_2015_PathoCNVs Wapner_2012_PathoCNVs DDD_2016_PathoCNVs ACMG_2013_PathoCNVs ; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_binomial.TBRden_binomial_annotation_test.results.txt | \
      cut -f6
    done
    for anno in Petrovski2013_RVIS_1 Clingen2015_Haploinsufficient Wang2015_essential Gencode_proteinCoding; do
      tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/unfiltered_annotations/${anno}/${group}_ANY_CNV_ANY_FILTER.${anno}.unfiltered_t.TBRden_ttest_annotation_test.results.txt | \
      cut -f6
    done
  done | paste -s
done >> ${WRKDIR}/data/plot_data/signif_loci_positive_control_enrichments.upperCI.txt

#####Prepare annotation files for functional enrichment analyses
#Make master bin file
zcat ${WRKDIR}/analysis/BIN_CNV_pileups/DD.DEL.TBRden_binned_pileup.bed.gz | \
cut -f1-4 | sed '1d' > ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.master_bins.bed
#Subset master bin file to autosomes only
zcat ${WRKDIR}/data/annotations/GRCh37.master_bins.bed.gz | fgrep -v "X" | fgrep -v "Y" | sed '1d' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
gzip -f ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed
#Exons (all) count
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_count.bed
#Exons (all) nonredundant bases
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.merged.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exonic_bases.bed
#Exons (all) nearest
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL\|^X\|^Y\|^M" ${SFARI_ANNO}/gencode/gencode.v25lift37.exons_and_UTRs.bed | \
  sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | awk -v OFS="\t" '{ if ($5==".") $5=="NA"; print }' > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exons_nearest.bed
#Exons (protein-coding) count
bedtools intersect -c -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed ) | \
sort -Vk1,1 -k2,2n -k3,3n > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_proteinCoding_count.bed
#Exons (protein-coding) nonredundant bases
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.exonic_proteinCoding_bases.bed
#Exons (protein-coding) nearest
bedtools closest -d -t first -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-b <( grep -v -e "^GL\|^X\|^Y\|^M" ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed | \
  sort -Vk1,1 -k2,2n -k3,3n ) | \
sort -Vk1,1 -k2,2n -k3,3n | cut -f1-4,9 > ${WRKDIR}/data/annotations/GRCh37.autosomes.exons_proteinCoding_nearest.bed
#GC content
bedtools nuc -fi ${h37} -bed ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz | \
cut -f1-5 | sed '1d' > ${WRKDIR}/data/annotations/GRCh37.autosomes.GC_pct.bed
#ENCODE DNAseH1 - Any
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_raw.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_raw_coverage.bed
#ENCODE DNAseH1 - 50% Cell Types
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_50pct_celltypes.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_50pct_coverage.bed
#ENCODE DNAseH1 - 90% Cell Types
bedtools coverage -b ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz \
-a <( grep -v -e "^GL" ${SFARI_ANNO}/noncoding/ENCODE2012_DNAse1_90pct_celltypes.bed ) | \
awk -v OFS="\t" '{ print $1, $2, $3, $4, $NF }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.DNAseH1_90pct_coverage.bed
#GM12878 histone marks
for mark in H3k27ac H3k4me1 H3k4me2 H3k4me3 H2az H3k36me3 H3k9ac H3k9me3; do
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/noncoding/ENCODE/wgEncodeBroadHistoneGm12878${mark}StdSig.bigWig \
  <( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab
  paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab ) > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.bed
  rm ${WRKDIR}/data/annotations/GRCh37.autosomes.${mark}.tab
done
#ENCODE 100mer alignability & 35bp uniqueness
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeCrgMapabilityAlign100mer.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Align_100mer.tab
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeDukeMapabilityUniqueness35bp.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Unique_35bp.tab
#GM12878 nucleosome map
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/wgEncodeSydhNsomeGm12878Sig.bigWig \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Nucleosome_Density.tab
#PhyloP 100-way conservation
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
/scratch/miket/rlc47temp/DOSAGE/hg19.100way.phyloP100way.bw \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.PhyloP_Conservation.tab
#Replication timing
bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o mean -c 4 \
-b <( fgrep -v "X" /scratch/miket/rlc47temp/DOSAGE/hg19_repTiming.bed | fgrep -v "Y" ) | \
awk -v OFS="\t" '{ if ($4==".") $4="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
${WRKDIR}/data/annotations/GRCh37.autosomes.RepTiming.bed
#Recombination frequency (deCODE)
${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
${SFARI_ANNO}/noncoding/ENCODE/SexAveraged.bw \
<( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab
paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
<( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab ) > \
${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.bed
rm ${WRKDIR}/data/annotations/GRCh37.autosomes.Recomb_Freq.tab
#GTEx Expression
#Curate
zcat ${SFARI_ANNO}/misc/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | \
sed -n '3p' | sed 's/\t/\n/g' | sed -e 's/\-//g' -e 's/(//g' -e 's/)//g' | sed 's/ \+/_/g' | \
sed '1,2d' | paste -s | paste <( echo -e "chr\tstart\tend\tgene" ) - > \
${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed
zcat ${SFARI_ANNO}/misc/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | \
sed '1,3d' | cut -f2- | sort -k1,1 | join -t $'\t' -1 1 -2 1 \
<( awk -v OFS="\t" '{ print $4, $1, $2, $3 }' ${SFARI_ANNO}/gencode/gencode.GRCh37.gene_boundaries.bed | sort -k1,1 ) - | \
awk -v OFS="\t" '{ if ($1!="7SK" && $1!="5S_rRNA" && $1!="Y_RNA") print $2, $3, $4, $1, $0 }' | \
cut -f 1-4,9- | sort -Vk1,1 -k2,2n -k3,3n >> \
${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed
while read tissue; do #master matrix of all genes
  echo ${tissue}
  col=$( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | fgrep -w ${tissue} | cut -f2 )
  cut -f1-4,${col} ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  sed '1d' > ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed
  gzip -f ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
while read tissue; do #Genes with FPKM > 50
  echo ${tissue}
  zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | \
  awk '{ if ($5>=50) print $0 }' > \
  ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed
  gzip -f ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.highExpressorGenes.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
while read tissue; do #curate binwise expression levels
  echo ${tissue}
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o mean -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_avgExpression_${tissue}.bed
  bedtools map -a ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed.gz -o max -c 5 \
  -b <( zcat ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.${tissue}.bed.gz | grep -v -e '^X\|^Y\|^M' ) | \
  awk -v OFS="\t" '{ if ($5==".") $5="NA"; print }' | sort -Vk1,1 -k2,2n -k3,3n > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.GTEx_maxExpression_${tissue}.bed
done < <( head -n1 ${SFARI_ANNO}/misc/GTEx_expression/GTEx_Expression.master_matrix.cleaned.bed | \
  cut -f4- | sed 's/\t/\n/g' | sed '1d' )
#Schmitt compartment - Fibroblast & LCL
for cell in BL AO AD tro SX SB RV PO PA OV npc msc mes LV LI LG imr90 HC h1 GM12878 CO; do
  ${SFARI_ANNO}/noncoding/ENCODE/bigWigAverageOverBed \
  ${SFARI_ANNO}/TADs/Schmitt2016/Compartment_primary_cohort/${cell}.pc.bw \
  <( awk -v OFS="\t" '{ print "chr"$0, NR }' ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed ) \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab
  paste ${WRKDIR}/data/annotations/GRCh37.autosomes.master_bins.bed \
  <( awk '{ print $5 }' ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab ) > \
  ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.bed
  rm ${WRKDIR}/data/annotations/GRCh37.autosomes.${cell}_compartment.tab
done

#####Get noncoding TBR fold-enrichments for noncoding CNVs
for group in ANY_DISEASE DD_SCZ CNCR; do
  echo ${group}
  for CNV in CNV DEL DUP; do
    tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/TBRs_MERGED/${group}_${CNV}_noncoding.TBRs_MERGED.filtered_binomial.TBRden_binomial_annotation_test.results.txt | \
    awk '{ print $2/$3 }'
    tail -n1 ${WRKDIR}/analysis/hotspot_enrichments/filtered_annotations/TBRs_MERGED/${group}_${CNV}_noncoding.TBRs_MERGED.filtered_binomial.TBRden_binomial_annotation_test.results.txt | \
    cut -f5,6
  done | paste -s 
done | paste - - > ${WRKDIR}/data/plot_data/Noncoding_rCNV_hotspot.noncoding_TBR_enrichments.txt

#####Run protein-coding exon pileups for noncoding rCNVs in flanks
#Prep list of 40kb flanks
distance=40000
awk -v OFS="\t" -v d=${distance} '{ print $1, $2-d, $3+d, $4 }' \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed | \
awk -v OFS="\t" '{ if ($2<0) $2=0; if ($3<1) $3=1; print }' > \
${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.40kb_windows.bed
#Launch tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/ ]; then
      rm -rf ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/
    fi
    mkdir ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_exon_flanking_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -n -z -t upper \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    -p ${group}_${CNV}_noncoding_exonFlanks \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.40kb_windows.bed \
    ${WRKDIR}/analysis/EXON_CNV_burdens/${group}_${CNV}_noncoding_exonFlanks/"
  done
done

#####Run TBRden TBR CNV pileups
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/ ]; then
      rm -rf ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/
    fi
    mkdir ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBR_pileup_all \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_pileup.sh -z \
    -o ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_all_TBR_pileups.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBR_pileup_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_pileup.sh -z \
    -o ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_noncoding_TBR_pileups.bed \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz"
  done
done

#####Run TBRden TBR CNV burden tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    if [ -e ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ ]; then
      rm -rf ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/
    fi
    mkdir ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/
    case ${group} in
      DD)
        color="green"
        ;;
      SCZ)
        color="purple"
        ;;
      DD_SCZ)
        color="blue"
        ;;
      CNCR)
        color="orange"
        ;;
    esac
    #Parallelize (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_all \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/TBR_CNV_pileups/CTRL_${CNV}_TBR_pileups/CTRL_${CNV}_all_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_all_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ \
    ${group}_${CNV}_all 0.000009072764 ${color}"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_noncoding \
    "${WRKDIR}/bin/rCNVmap/bin/TBRden_test.R \
    ${WRKDIR}/analysis/TBR_CNV_pileups/CTRL_${CNV}_TBR_pileups/CTRL_${CNV}_noncoding_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_pileups/${group}_${CNV}_TBR_pileups/${group}_${CNV}_noncoding_TBR_pileups.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/ \
    ${group}_${CNV}_noncoding 0.000009072764 ${color}"
  done
done

#####Run direct permutation tests of TBRs 
#Launch tests
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    #Parallelize intersections (LSF)
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_TBR_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -n -z -t upper \
    -e ${SFARI_ANNO}/gencode/gencode.v25lift37.protein_coding_exons.no_ASmerged.bed \
    -p ${group}_${CNV}_TBR_burdens_noncoding \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/"
    bsub -q short -sla miket_sc -u nobody -J ${group}_${CNV}_TBRden_TBR_burden \
    "${WRKDIR}/bin/rCNVmap/bin/direct_burden_test.sh -d 5 -N 1000 -z -t upper \
    -p ${group}_${CNV}_TBR_burdens_all \
    ${WRKDIR}/data/CNV/CNV_MASTER/CTRL.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz \
    ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz \
    ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/"
  done
done

#####Aggregate results of TBR noncoding rCNV burden tests
#Print header
for dummy in 1; do
  echo -e "#chr\tstart\tend\tTBR_ID\ttissues"
  for group in DD SCZ DD_SCZ CNCR; do
    for CNV in DEL DUP CNV; do
      echo -e "${group}_${CNV}_p_obs\n${group}_${CNV}_p_perm"
    done
  done | paste -s
done | paste -s > ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
#Gather p-values
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    zcat ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/${group}_${CNV}_noncoding.TBRden_results.bed.gz | \
    sed '1d' | awk '{ print $NF }' | paste -s
    zcat ${WRKDIR}/analysis/TBR_CNV_burdens/${group}_${CNV}_TBR_burdens/${group}_${CNV}_TBR_burdens_noncoding.TBRden_direct_test_results.bed.gz | \
    sed '1d' | awk '{ print $(NF-1) }' | paste -s
  done
done > ${TMPDIR}/matrix_pretranspose.txt
#Transpose matrix
Rscript -e "write.table(t(read.table(\"${TMPDIR}/matrix_pretranspose.txt\",header=F)),\"${TMPDIR}/matrix_posttranspose.txt\",col.names=F,row.names=F,sep=\"\\t\")"
paste <( zcat ${WRKDIR}/data/unfiltered_annotations/TBRs_MERGED.boundaries.bed.gz | fgrep -v "#" | cut -f1-5 ) \
${TMPDIR}/matrix_posttranspose.txt >> ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
gzip -f ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt
#Copy to plot directory
cp ${WRKDIR}/analysis/TBR_CNV_burdens/MASTER_TBR_rCNV_burdens.results.txt.gz ${WRKDIR}/data/plot_data/



