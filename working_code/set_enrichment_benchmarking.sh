#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Master annotation curation code

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

#####Simulate test intervals
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals
#Copy genome and restrict to autosomes
grep -e '^[1-9]' /data/talkowski/rlc47/src/GRCh37.genome > \
${WRKDIR}/data/misc/GRCh37_autosomes.genome
#Simulate 100 sets of intervals with various set of parameters
for n in 10 100 1000 10000 100000; do
  for i in $( seq -w 001 100 ); do
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
    done < <( echo -e "50\t10\n500\t100\n5000\t1000\n50000\t10000" )
  done
done

#####Launch simple annotation shuffling permutation test (10k permutations each)
#Test rCNV sets: all control deletions vs all germline deletions
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_rCNV
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_urCNV
for n in 10 100 1000 10000 100000; do
  for i in $( seq -w 001 100 ); do
    #size = 50bp, stdev = 10bp, n=10, 100, 1000, 10000, 100000
    #size = 500bp, stdev = 100bp, n=10, 100, 1000, 10000, 100000
    #size = 5000bp, stdev = 1000bp, n=10, 100, 1000, 10000, 100000
    #size = 50000bp, stdev = 10000bp, n=10, 100, 1000, 10000, 100000
    while read size sd; do
      bsub -q short -sla miket_sc -u nobody -J permutation_test_rCNV_${size}bp_${sd}bp_x${n}_i${i} \
      "${WRKDIR}/bin/rCNVmap/bin/annoSet_permutation_test.sh -N 10000 \
       -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
       -o ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_rCNV/results_${size}bp_${sd}bp_x${n}_i${i}.txt \
       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.GRCh37.all.bed.gz \
       ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.DEL.GRCh37.all.bed.gz \
       ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed.gz \
       ${WRKDIR}/data/misc/GRCh37_autosomes.genome"
      bsub -q short -sla miket_sc -u nobody -J permutation_test_urCNV_${size}bp_${sd}bp_x${n}_i${i} \
      "${WRKDIR}/bin/rCNVmap/bin/annoSet_permutation_test.sh -N 10000 \
       -x /data/talkowski/rlc47/src/GRCh37_Nmask.bed \
       -o ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_urCNV/results_${size}bp_${sd}bp_x${n}_i${i}.txt \
       ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.DEL.urCNVs.GRCh37.all.bed.gz \
       ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.DEL.urCNVs.GRCh37.all.bed.gz \
       ${WRKDIR}/analysis/benchmarking/set_enrichments/simulated_intervals/intervals_${size}bp_${sd}bp_x${n}_i${i}.bed.gz \
       ${WRKDIR}/data/misc/GRCh37_autosomes.genome"
    done < <( echo -e "50\t10\n500\t100\n5000\t1000\n50000\t10000" )
  done
done
