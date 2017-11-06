# for CNV in DEL DUP; do
#   for pheno in GERM NEURO NDD PSYCH SOMA; do  
#     cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_${VF}_${context}.geneScore_Bonferroni_sig.FINAL.genes.list
#   done | sort | u
# done > ${TMPDIR}/germ.list

# for pheno in CNCR; do
#   for CNV in DEL DUP; do
#     cat /data/talkowski/Samples/rCNVmap/analysis/perGene_burden/signif_genes/merged/${pheno}_${CNV}_E4_exonic.geneScore_Bonferroni_sig.union.genes.list
#   done
# done | sort | uniq | fgrep -wf ${TMPDIR}/germ.list

#Germline risk estimates
#Large loci
for CNV in DEL DUP; do
  for wrapper in 1; do
    bedtools intersect -wb -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f7 | sort | uniq | wc -l
    bedtools intersect -wb -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f7 | sort | uniq | wc -l
  done | paste -s | awk -v OFS="\t" '{ print $1/63629, $2/38628 }' | awk -v OFS="\t" '{ print $1, $2, $1-$2 }'
done
#Genes
for CNV in DEL DUP; do
  for wrapper in 1; do
    fgrep -wf ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/MasterPhenoGroups_${CNV}_E4_${context}.geneScore_Bonferroni_sig.FINAL.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed |
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f4 | sort | uniq | wc -l
    fgrep -wf ${TMPDIR}/GERM_all.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed |
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f4 | sort | uniq | wc -l
  done | paste -s | awk -v OFS="\t" '{ print $1/63629, $2/38628 }' | awk -v OFS="\t" '{ print $1, $2, $1-$2 }'
done
#Noncoding
for CNV in DEL DUP; do
  if [ ${CNV} == "DEL" ]; then
    filt=haplosufficient
  else
    filt=noncoding
  fi
  for wrapper in 1; do
    zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    cut -f4 | sort | uniq | wc -l
    zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    cut -f4 | sort | uniq | wc -l
  done | paste -s | awk -v OFS="\t" '{ print $1/63629, $2/38628 }' | awk -v OFS="\t" '{ print $1, $2, $1-$2 }'
done

#Germline relative risk estimates per estimate (unweighted for CNVs overlapping multiple elements)
#Large loci
for CNV in DEL DUP; do
  for wrapper in 1; do
    paste <( bedtools intersect -c -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f4 ) \
    <( bedtools intersect -c -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f4 ) | awk -v OFS="\t" '{ print $1, $2, 63629-$1, 388628-$2 }' | \
    awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'
  done
done | sort -nrk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#Genes
for CNV in DEL DUP; do
  for wrapper in 1; do
    paste <( fgrep -wf ${TMPDIR}/GERM_all.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed |
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f5 ) \
    <( fgrep -wf ${TMPDIR}/GERM_all.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed |
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f5 ) | awk -v OFS="\t" '{ print $1+0.0001, $2+0.0001, 63629-$1, 388628-$2 }' | \
    awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'
  done 
done | sort -nrk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'
#Noncoding
for CNV in DEL DUP; do
  if [ ${CNV} == "DEL" ]; then
    filt=haplosufficient
  else
    filt=noncoding
  fi
  for wrapper in 1; do
    paste <( zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    awk '{ print $NF }' ) \
    <( zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -c -a - \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    awk '{ print $NF }' ) | awk -v OFS="\t" '{ print $1, $2, 63629-$1, 388628-$2 }' | \
    awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'
  done
done | sort -nrk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'


#Germline relative risk estimates per estimate in total for each category
#Large loci
for CNV in DEL DUP; do
  for wrapper in 1; do
    bedtools intersect -wb -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f7 | sort | uniq | wc -l
    bedtools intersect -wb -f 0.8 \
    -a ${WRKDIR}/analysis/large_CNV_segments/GERM/GERM_${CNV}_E2_all.signif.bed \
    -b ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E2.GRCh37.all.bed.gz | \
    cut -f7 | sort | uniq | wc -l
  done | paste -s
done | awk -v OFS="\t" '{ case_sum+=$1; ctrl_sum+=$2}END{ print case_sum, ctrl_sum }' | \
awk -v OFS="\t" '{ print $1, $2, 63629-$1, 388628-$2 }' | awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'
#Genes
for CNV in DEL DUP; do
  for wrapper in 1; do
    fgrep -wf ${TMPDIR}/GERM_all.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed |
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f4 | sort | uniq | wc -l
    fgrep -wf ${TMPDIR}/GERM_all.genes.list \
    ${WRKDIR}/data/master_annotations/gencode/gencode.v19.exons.protein_coding.bed |
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.coding.bed.gz | \
    cut -f4 | sort | uniq | wc -l
  done | paste -s
done | awk -v OFS="\t" '{ case_sum+=$1; ctrl_sum+=$2}END{ print case_sum, ctrl_sum }' | \
awk -v OFS="\t" '{ print $1, $2, 63629-$1, 388628-$2 }' | awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'
#Noncoding
for CNV in DEL DUP; do
  if [ ${CNV} == "DEL" ]; then
    filt=haplosufficient
  else
    filt=noncoding
  fi
  for wrapper in 1; do
    zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/GERM/GERM.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    cut -f4 | sort | uniq | wc -l
    zcat ${WRKDIR}/analysis/perAnno_burden/cleaned_noncoding_loci.wAnnotations.forManualCuration.txt.gz | \
    cut -f1-3,5-9 | fgrep ${CNV} | \
    bedtools intersect -u -b - \
    -a ${WRKDIR}/data/CNV/CNV_MASTER/CTRL/CTRL.${CNV}.E4.GRCh37.${filt}.bed.gz | \
    cut -f4 | sort | uniq | wc -l
  done | paste -s | awk -v OFS="\t" '{ print $1/63629, $2/38628 }' | awk -v OFS="\t" '{ print $1, $2, $1-$2 }'
done | awk -v OFS="\t" '{ case_sum+=$1; ctrl_sum+=$2}END{ print case_sum, ctrl_sum }' | \
awk -v OFS="\t" '{ print $1, $2, 63629-$1, 388628-$2 }' | awk '{ print ($1/($1+$2)) / ($3/($3+$4)) }'









