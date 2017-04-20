#Trial code for new 1kb rCNV hotspot analysis

#####Set parameters
export WRKDIR=/data/talkowski/Samples/rCNVmap
source ${WRKDIR}/bin/rCNVmap/misc/rCNV_code_parameters.sh

#####Segment genome into 1kb bins
while read chr length; do
  paste <( seq 0 1000 $((${length}-1000)) ) \
  <( seq 1000 1000 ${length} ) | \
  awk -v OFS="\t" -v chr=${chr} '{ print chr, $0 }'
done < /data/talkowski/rlc47/src/GRCh37.genome > \
${WRKDIR}/data/misc/GRCh37.1kb_bins.all.bed
gzip -f ${WRKDIR}/data/misc/GRCh37.1kb_bins.all.bed

#####Filter 1kb bins
#Restrict to autosomes
zcat ${WRKDIR}/data/misc/GRCh37.1kb_bins.all.bed.gz | grep -ve '^X\|^Y\|^M\|^GL' > \
${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.bed
gzip -f ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.bed
#Remove N-masked regions
bedtools intersect -v -a ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.bed.gz \
-b /data/talkowski/rlc47/src/GRCh37_Nmask.bed > \
${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.bed
gzip -f ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.bed
#Remove bins within 1Mb of centromere or telomere
awk -v OFS="\t" '{ print $1, 0, 1000000"\n"$1, $2-1000000, $2 }' \
/data/talkowski/rlc47/src/GRCh37.genome | \
awk -v OFS="\t" '{ if ($2<0) $2=0; if ($3<1) $3=1; print }' | \
bedtools intersect -v -b - -a ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.bed.gz | \
awk -v OFS="\t" '{ print $0, "Bin_"$1"_"NR }' > \
${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.Blacklisted.bed
gzip -f ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.Blacklisted.bed

#####Run TBRden CNV pileups
for group in CTRL DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    echo -e "${group} ${CNV} all"
    bedtools intersect -c -wa \
    -a ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.Blacklisted.bed.gz \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.bed.gz > \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.GRCh37_1kb_bins.pileup_all.bed
    echo -e "${group} ${CNV} coding"
    bedtools intersect -c -wa \
    -a ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.Blacklisted.bed.gz \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.coding.bed.gz > \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.GRCh37_1kb_bins.pileup_coding.bed
    echo -e "${group} ${CNV} noncoding"
    bedtools intersect -c -wa \
    -a ${WRKDIR}/data/misc/GRCh37.1kb_bins.autosomes.Nmasked.Blacklisted.bed.gz \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/${group}.${CNV}.GRCh37.noncoding.bed.gz > \
    ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.GRCh37_1kb_bins.pileup_noncoding.bed
  done
done
gzip -f ${WRKDIR}/analysis/BIN_CNV_pileups/*.GRCh37_1kb_bins.pileup_*.bed

#####Run preliminary genome-wide Fisher's exact tests versus case & control medians
for group in DD SCZ DD_SCZ CNCR; do
  # if [ -e ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL ]; then
  #   rm -rf ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL
  # fi
  for CNV in DEL DUP CNV; do
    for filt in all coding noncoding; do
      bsub -q big -R 'rusage[mem=16000]' -M 16000 -v 20000 -sla miket_sc -u nobody -J ${group}_${CNV}_FisherGenome \
      "${WRKDIR}/bin/rCNVmap/bin/FisherGenome.R -z \
      ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.GRCh37_1kb_bins.pileup_${filt}.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_pileups/CTRL.${CNV}.GRCh37_1kb_bins.pileup_${filt}.bed.gz \
      ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.GRCh37_1kb_bins.FisherGenome.bed"
    done
  done
done

#####Count number of bins examined, tested, Fisher significant, and permuted significant
for group in DD SCZ DD_SCZ CNCR; do
  for CNV in DEL DUP CNV; do
    for filt in all coding noncoding; do
      for dummy in 1; do
        echo -e "${group}\t${CNV}\t${filt}"
        #Bins examined
        zcat ${WRKDIR}/analysis/BIN_CNV_pileups/${group}.${CNV}.GRCh37_1kb_bins.pileup_${filt}.bed.gz | wc -l
        #Bins Fisher tested
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.GRCh37_1kb_bins.FisherGenome.bed.gz | \
        fgrep -v "#" | awk '{ if ($7!="NA") print $0 }' | wc -l
        #Bins Fisher significant
        zcat ${WRKDIR}/analysis/BIN_CNV_burdens/${group}_vs_CTRL/${group}_vs_CTRL_${CNV}_${filt}.GRCh37_1kb_bins.FisherGenome.bed.gz | \
        fgrep -v "#" | awk '{ if ($7!="NA" && $7<=0.05) print $0 }' | wc -l
      done | paste -s
    done
  done
done


