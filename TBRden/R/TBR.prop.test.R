#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#TBR.prop.test: function to perform TBR proportion test

TBR.prop.test <- function(TBR.CTRL,  #Size-normalized CNV score from TBRden_pileup.sh for controls in TBR
                          TAD.CTRL,  #Size-normalized CNV score from TBRden_pileup.sh for controls in TAD
                          TBR.CASE,  #Size-normalized CNV score from TBRden_pileup.sh for cases in TBR
                          TAD.CASE   #Size-normalized CNV score from TBRden_pileup.sh for cases in TAD
){
  #Coerce all input into vector of numerics
  scores <- as.numeric(c(TBR.CTRL,TAD.CTRL,TBR.CASE,TAD.CASE))

  #Compute control proportion
  if(!(any(is.na(c(TBR.CTRL,TAD.CTRL))))){
    if(TBR.CTRL==0 & TAD.CTRL==0){
      prop.CTRL=NA
    }else{
      prop.CTRL <- scores[1]/sum(scores[1:2])
    }
  }else{
    prop.CTRL=NA
  }

  #Compute control proportion
  if(!(any(is.na(c(TBR.CASE,TAD.CASE))))){
    if(TBR.CASE==0 & TAD.CASE==0){
      prop.CASE=NA
    }else{
      prop.CASE <- scores[3]/sum(scores[3:4])
    }
  }else{
    prop.CASE=NA
  }

  #Unable to run test if both TBR and TAD scores are 0 for either group
  if(any(is.na(c(prop.CTRL,prop.CASE)))){
    D <- NA
  }else{
    D <- prop.CASE - prop.CTRL
  }

  #Return proportions from case & control and difference
  return(c(prop.CTRL,prop.CASE,D))

}
