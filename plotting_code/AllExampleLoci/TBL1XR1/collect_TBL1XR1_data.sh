#NLGN1-TBL1XR1 compartment boundary
#GERM DEL: INT, SKEL, MSC
#CNCR DEL: CMSK
#GERM DUP: ASD, SEIZ, GRO
#CNCR DUP: CSKN, CMSK, CBST, CLNG
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL DUP; do
      bedtools intersect -wb -a <( echo -e "3\t175748655\t176675143" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )


#TBL1XR1 (distal)
#GERM: nothing
#CNCR DUP: CSKN, CMSK, CGEN, CLNG
#CNCR DEL: CMSK
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL DUP; do
      bedtools intersect -wb -a <( echo -e "3\t176930083\t177638939" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )


#NLGN1 (proximal)
#GERM DEL: NDD/DD with a smattering of other phenos
#CNCR DEL: CGST, CMSK
#GERM DUP: ASD (over NAALADL2)
#CNCR DUP: CSKN, CBST, CLNG, CGEN, CHNK
while read pheno nsamp; do
  for dummy in 1; do
    echo -e "${pheno}\t${nsamp}"
    for CNV in DEL DUP; do
      bedtools intersect -wb -a <( echo -e "3\t174069784\t175189032" ) \
      -b ${WRKDIR}/data/CNV/CNV_MASTER/${pheno}/${pheno}.${CNV}.E4.GRCh37.haplosufficient.bed.gz | wc -l
    done
  done | paste -s
done < <( fgrep -v "#" ${WRKDIR}/bin/rCNVmap/misc/analysis_group_HPO_mappings.list | \
          awk -v OFS="\t" '{ if ($2!="CTRL") print $1, $NF }' )


#Compartments:
#Cortex & hippocampus proximal (173-177): moderate B, strongest B is over deletion hotspot
#Cortex & hippocampus distal (177-182): moderate A, gets much stronger in 1-2Mb beyond distal region
#IMR90, Psoas muscle, lung proximal: weak/moderate B
#IMR90, Psoas muscle, lung distal: extremely strong A for just 1Mb (~right over TBL1XR1), then reverts back to strong B




#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Reinitialize directory if exists
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/ ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/
fi
if ! [ -e ${WRKDIR}/data/plot_data/ExampleLocusPlots/TBL1XR1 ]; then
  mkdir ${WRKDIR}/data/plot_data/ExampleLocusPlots/TBL1XR1
fi

#####Gather compartment data in ~5kb bins
#Write 5kb bins
paste <( seq 170300000 5000 180695000 ) <( seq 170305000 5000 180700000 ) | \
awk -v OFS="\t" '{ print "chr3", $1, $2, "chr3:"$1"-"$2 }' > \
${WRKDIR}/data/plot_data/ExampleLocusPlots/TBL1XR1/TBL1XR1_locus.hg19.5kb_bins.bed
#Iterate over tissues and compute compartment scores
while read code tissue; do
  echo ${tissue}
  ${WRKDIR}/bin/utils/bigWigAverageOverBed \
  -bedOut=${WRKDIR}/data/plot_data/ExampleLocusPlots/TBL1XR1/${tissue}.compartments_5kbBins.bed \
  /data/talkowski/Samples/SFARI/ASC_analysis/annotations/TADs/Schmitt2016/Compartment_primary_cohort/${code}.pc.bw \
  ${WRKDIR}/data/plot_data/ExampleLocusPlots/TBL1XR1/TBL1XR1_locus.hg19.5kb_bins.bed \
  /dev/null
done < <( echo -e "CO\tcortex\nHC\thippocampus\nimr90\tIMR90\nLG\tlung\nPO\tmuscle" )









