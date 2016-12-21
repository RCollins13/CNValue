#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#TBR.fischer.test: function to perform TBR Fisher's exact test

TBR.fisher.test <- function(TBR.CTRL,  #Count of CNVs overlapping TBR in controls
                            TAD.CTRL,  #Count of CNVs overlapping TAD in controls
                            TBR.CASE,  #Size-normalized CNV score from TBRden_pileup.sh for cases in TBR
                            TAD.CASE   #Size-normalized CNV score from TBRden_pileup.sh for cases in TAD
){
  #Coerce all input into vector of numerics
  scores <- as.numeric(c(TBR.CTRL,TAD.CTRL,TBR.CASE,TAD.CASE))

  #Compute odds ratio
  if(!(any(is.na(scores)))){
    if(TAD.CTRL==0 | TAD.CASE==0){
      OR <- NA
    }else{
      OR <- (scores[3]/scores[4])/(scores[1]/scores[2])
    }
  }else{
    OR <- NA
  }

  #Run Fisher's exact test
  if(!(any(is.na(scores)))){
    p <- fisher.test(matrix(scores[c(3:4,1:2)],nrow=2),
                     alternative="greater")$p.value
  }else{
    p <- NA
  }

  #Return proportions from case & control and difference
  return(c(p,OR))
}
