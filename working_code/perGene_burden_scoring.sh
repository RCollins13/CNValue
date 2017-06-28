#!/usr/bin/env bash

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Code to run all rCNV burden scoring per gene

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if [ -e ${WRKDIR}/data/perGene_burden ]; then
  rm -rf ${WRKDIR}/data/perGene_burden
fi
mkdir ${WRKDIR}/data/perGene_burden

#####Create subdirectories
for pheno in GERM NEURO SOMA CNCR; do
  if [ -e ${WRKDIR}/data/perGene_burden/${pheno} ]; then
    rm -rf ${WRKDIR}/data/perGene_burden/${pheno}
  fi
  mkdir ${WRKDIR}/data/perGene_burden/${pheno}
done

#####Submit burden data collection for GERM, NEURO, and SOMA at E2-N1 and CNCR at E2
#Germline
for pheno in GERM NEURO SOMA; do
  for CNV in CNV DEL DUP; do
    for VF in E2 E3 E4 N1; do
      #Exonic
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_exonic \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_exonic.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
      #Wholegene
      bsub -q short -sla miket_sc -u nobody -J ${pheno}_${CNV}_${VF}_perGene_burden_dataCollection_wholegene \
      "${WRKDIR}/bin/rCNVmap/bin/gather_geneScore_data.sh -W \
      -H ${WRKDIR}/data/misc/exons_boundaries_dictionary/ \
      -U /data/talkowski/Samples/rCNVmap/data/master_annotations/genelists/Gencode_v19_protein_coding.genes.list \
      -o ${WRKDIR}/data/perGene_burden/${pheno}/${pheno}_${CNV}_${VF}_wholegene.geneScore_data.txt \
      ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.${VF}.GRCh37.${filt}.bed.gz \
      ${WRKDIR}/data/master_annotations/gencode/gencode.v19.annotation.gtf \
      ${h37}"
    done
  done
done

