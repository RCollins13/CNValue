#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code for benchmarking set enrichment tests

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

#####Prepare annotation directory tree
mkdir ${WRKDIR}/analysis/benchmarking
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments

######Simulate 1,000 sets of intervals with various set of parameters
#Copy genome and restrict to autosomes
grep -e '^[1-9]' /data/talkowski/rlc47/src/GRCh37.genome > \
${WRKDIR}/data/misc/GRCh37_autosomes.genome
#Create directory
if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals ]; then
  rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals
fi
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals
#Simulate test intervals
for n in 10 100 1000 10000 100000 1000000; do
  for i in $( seq -w 0001 1000 ); do
    #size = 50bp, stdev = 10bp, n=10, 100, 1000, 10000, 100000
    #size = 500bp, stdev = 100bp, n=10, 100, 1000, 10000, 100000
    #size = 5000bp, stdev = 1000bp, n=10, 100, 1000, 10000, 100000
    #size = 50000bp, stdev = 10000bp, n=10, 100, 1000, 10000, 100000
    while read size sd; do
      bsub -q short -sla miket_sc -u nobody -J generate_intervals_${size}bp_${sd}bp_x${n}_i${i} \
      "${WRKDIR}/bin/rCNVmap/bin/simulate_intervals.sh -z -s ${size} -d ${sd} -N ${n} \
       -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
       -o ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed \
       ${WRKDIR}/data/misc/GRCh37_autosomes.genome"
    done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
  done
done

#####Test rCNVs, vrCNVs, and sCNV for CNV/DEL/DUP for GERM and CNCR
#1k tests per simulated set of parameters
#10k permutations per test
for set in rCNV vrCNV sCNV; do
  if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}/${pheno}
      #size = 5bp, stdev = 1bp, n=10, 100, 1000, 10000, 100000, 1000000
      #size = 50bp, stdev = 10bp, n=10, 100, 1000, 10000, 100000, 1000000
      #size = 500bp, stdev = 100bp, n=10, 100, 1000, 10000, 100000, 1000000
      #size = 5000bp, stdev = 1000bp, n=10, 100, 1000, 10000, 100000, 1000000
      #size = 50000bp, stdev = 10000bp, n=10, 100, 1000, 10000, 100000, 1000000
      for n in 10 100 1000 10000 100000 1000000; do
        while read size sd; do
          #Code to launch simulations per all 1k test sets
            bsub -q short -sla miket_sc -u nobody -J ${set}_${CNV}_${pheno}_annoSet_permutation_test_${size}bp_${sd}bp_x${n} \
            "${WRKDIR}/bin/rCNVmap/analysis_scripts/run_set_enrichment_permutations.sh \
             ${set} ${CNV} ${pheno} ${n} ${size} ${sd}"
           done
        done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
      done
    done
  done
done

#####Collect results from annotation set shuffling permutation tests
# mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results
# for n in 10 100 1000 10000 100000; do
#   while read size sd; do
#     for set in rCNV urCNV; do
#       for i in $( seq -w 001 100 ); do
#         fgrep -v "#" ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/results_${size}bp_${sd}bp_x${n}_i${i}.txt | awk '{ print $NF }'
#        done > ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}_results_${size}bp_${sd}bp_x${n}.txt
#     done
#   done < <( echo -e "50\t10\n500\t100\n5000\t1000\n50000\t10000" )
# done



