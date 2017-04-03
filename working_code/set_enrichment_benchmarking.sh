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
mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments

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
        done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
      done
    done
  done
done

#####Collect results from annotation set shuffling permutation tests
if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results ]; then
  rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results
fi
mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results
for VF in E2 E3 E4 N1; do
  if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/set_enrichments/results/${VF}/${CNV}/${pheno}
      for n in 10 100 1000 10000 100000 1000000; do
        while read size sd; do
          bsub -q short -sla miket_sc -J collectPermutations_${pheno}_${CNV}_${VF}_${size}bp_x${n} -u nobody \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_set_enrichment_permutations.sh \
          ${VF} ${CNV} ${pheno} ${n} ${size} ${sd}"
        done < <( echo -e "5\t1\n50\t10\n500\t100\n5000\t1000\n50000\t10000" )
      done
    done
  done
done

######Simulate 1,000 random gene sets with various set of parameters
#Create directory
if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/simulated_sets ]; then
  rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/simulated_sets
fi
mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/simulated_sets
#Simulate sets
for n in 5 10 50 100 1000 5000 10000; do
  echo -e "\n\n\n${n}\n\n\n"
  mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/simulated_sets/genes_n${n}
  for i in $( seq -w 0001 1000 ); do
    echo ${i}
    shuf ${WRKDIR}/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list | head -n${n} > \
    ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/simulated_sets/genes_n${n}/geneSet_n${n}_${i}.list
  done
done

#####Generate filtered exon & boundary files for hard override (saves time)
mkdir ${WRKDIR}/data/misc/exons_boundaries_dictionary
#All exons
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
'{ if ($3=="exon") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
sed 's/^chr//g' | sed 's/\-/_/g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > \
${WRKDIR}/data/misc/exons_boundaries_dictionary/exons.bed
#All boundaries
fgrep -v "#" ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf | \
sed 's/gene_name/\t/g' | awk -v FS="\t" -v OFS="\t" \
'{ if ($3=="gene") print $1, $4, $5, $10 }' | sed 's/\;/\t/g' | \
awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' | tr -d "\"" | \
sed 's/^chr//g' | sed 's/\-/_/g' | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > \
${WRKDIR}/data/misc/exons_boundaries_dictionary/boundaries.bed

#####Test all VF filters for CNV/DEL/DUP for GERM and CNCR
#100 independent tests per gene set size
#1k permutations per test
for VF in E2 E3 E4 N1; do
  if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/permutation_testing_${VF}/${CNV}/${pheno}
      for n in 5 10 50 100 1000 5000 10000; do
        for W in 0 1; do
          #Code to launch simulations per all 1k test sets
          bsub -q short -sla miket_sc -u nobody -J ${VF}_${CNV}_${pheno}_geneSet_permutation_test_n${n}_W${W} \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/run_geneSet_enrichment_permutations.sh \
           ${VF} ${CNV} ${pheno} ${n} ${W}"
         done
      done
    done
  done
done

#####Collect results from gene set shuffling permutation tests
if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results ]; then
  rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results
fi
mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results
for VF in E2 E3 E4 N1; do
  if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF} ]; then
    rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}
  fi
  mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}
  for CNV in CNV DEL DUP; do
    if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV} ]; then
      rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}
    fi
    mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}
    for pheno in GERM CNCR; do
      if [ -e ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}/${pheno} ]; then
        rm -rf ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}/${pheno}
      fi
      mkdir ${WRKDIR}/analysis/benchmarking/geneSet_enrichments/results/${VF}/${CNV}/${pheno}
      for n in 5 10 50 100 1000 5000 10000; do
        for W in 0 1; do
          bsub -q short -sla miket_sc -J collectPermutations_${pheno}_${CNV}_${VF}_${size}bp_x${n} -u nobody \
          "${WRKDIR}/bin/rCNVmap/analysis_scripts/collect_set_enrichment_permutations.sh \
          ${VF} ${CNV} ${pheno} ${n} ${size} ${sd}"
        done
      done
    done
  done
done

















