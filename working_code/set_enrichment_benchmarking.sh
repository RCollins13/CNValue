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
  while read size sd; do
    bsub -q short -sla miket_sc -u nobody -J generate_intervals_${size}bp_${sd}bp_x${n} \
    "${WRKDIR}/bin/rCNVmap/analysis_scripts/generate_simulated_intervals.sh ${size} ${sd} ${n}"
  done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
done

#####Test all VF filters for CNV/DEL/DUP for GERM and CNCR
#1k independent tests per simulated set of parameters
#1k permutations per test
for VF in E2 E3 E4 N1; do
  if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${VF}/${CNV}/${pheno}
      for n in 10 100 1000 10000 100000 1000000; do
        while read size sd; do
          #Code to launch simulations per all 1k test sets
            bsub -q short -sla miket_sc -u nobody -J ${VF}_${CNV}_${pheno}_annoSet_permutation_test_${size}bp_${sd}bp_x${n} \
            "${WRKDIR}/bin/rCNVmap/analysis_scripts/run_set_enrichment_permutations.sh \
             ${VF} ${CNV} ${pheno} ${n} ${size} ${sd}"
           done
        done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
      done
    done
  done
done

#####Collect results from annotation set shuffling permutation tests
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results
for set in rCNV vrCNV sCNV; do
  if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}/${pheno}
      for n in 10 100 1000 10000 100000 1000000; do
        while read size sd; do
          for i in $( seq -w 0001 1000 ); do
            if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt ]; then
              fgrep -v "#" ${WRKDIR}/analysis/benchmarking/set_enrichments/permutation_testing_${set}/${CNV}/${pheno}/results_${size}bp_${sd}bp_x${n}_i${i}.txt | awk '{ print $NF }'
            fi
           done > ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${set}/${CNV}/${pheno}/${set}_results_${size}bp_${sd}bp_x${n}.txt
        done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
      done
    done
  done
done



