#!/usr/bin/env R

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#fitGeneModel: corrects observed dCNVs per gene by
# linear regression versus gene length and exonic base

fitGeneModel <- function(path,
                         outfile,
                         plot=T){
  #Sanity check input file
  if(!(file.exists(path))){
    stop(paste("Input file ",path," does not exist.",sep=""))
  }

  #Read data
  dat <- read.table(path,header=F)
  names(dat) <- c("gene","length","exon.bp","dCNV")

  #Fit multivariate linear model
  model <- lm(dCNV ~ length + exon.bp, data=dat)

  predict.lm(model,dat)





}
