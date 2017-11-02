#Get dup-sig genes
for pheno in GERM NDD NEURO; do
  cat ${WRKDIR}/analysis/perGene_burden/signif_genes/merged/${pheno}_DUP_E4_wholegene.geneScore_Bonferroni_sig.unique.genes.list
done | sort | uniq > ${TMPDIR}/candidate_genes.list

#Subset those which are haplosufficient
while read gene; do
  if [ $( fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/ExAC_haplosufficient.genes.list | sed 's/\-/\t/g' | cut -f1 | sed '/^$/d' | wc -l ) -gt 0 ]; then
    echo ${gene}
  fi
done < ${TMPDIR}/candidate_genes.list > \
${TMPDIR}/candidate_genes.haplosuff.list

#Subset to those which are highly expressed in at least one tissue
while read gene; do
  if [ $( fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/*Highly_Expressed.genes.list | sed 's/\-/\t/g' | cut -f1 | sed '/^$/d' | wc -l ) -gt 0 ]; then
    echo ${gene}
  fi
done < ${TMPDIR}/candidate_genes.haplosuff.list > \
${TMPDIR}/candidate_genes.haplosuff.high_expr.list

#Subset to those which are lowly expressed in at least one tissue
while read gene; do
  if [ $( fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/*Lowly_Expressed.genes.list | sed 's/\-/\t/g' | cut -f1 | sed '/^$/d' | wc -l ) -gt 0 ]; then
    echo ${gene}
  fi
done < ${TMPDIR}/candidate_genes.haplosuff.high_expr.list > \
${TMPDIR}/candidate_genes.haplosuff.high_expr.low_expr.list

#Iterate over those genes and print gene lists
while read gene; do
  echo -e "\n\n${gene}--"
  fgrep -w ${gene} ${WRKDIR}/data/master_annotations/genelists/* | sed 's/\-/\t/g' | cut -f1 | fgrep -w ${gene} | fgrep -v Expr | fgrep -v GTEx | fgrep -v Gencode | fgrep -v ExAC
done < ${TMPDIR}/candidate_genes.haplosuff.high_expr.low_expr.list

#Get missense z-scores for all genes
while read gene; do
  for wrapper in 1; do
    echo -e "${gene}"
    fgrep -w ${gene} ${WRKDIR}/data/misc/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt | awk '{ print $18 }'
  done | paste -s
done < ${TMPDIR}/candidate_genes.haplosuff.high_expr.low_expr.list


#Top candidates
# CD8A (best)
# CYTL1 (best)
# GEM (best)
# KCNJ12
# MMP17
# HGFAC


#Reasonable candidate genes with concerns about being subtelomeric or pericentromeric
# FAM110C
# MICALL2
# CHL1
# EMB
