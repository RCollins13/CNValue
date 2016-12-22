#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Add columns for cross-disease and cross-CNV comparisions from master merged burden testing

#Set params
options(scipen=1000)

#Read arguments
#first: input
#second: output
args <- commandArgs(trailingOnly=T)

#Read table
x <- read.table(args[1],comment.char="",header=T)

#Run
merged_results <- as.data.frame(t(sapply(1:nrow(x),function(i){
  vals <- x[i,-c(1:4)]
  names(vals) <- colnames(x)[-c(1:4)]
  #All CNVs, coding + noncoding
  MIN.CNV.all.obs_p <- min(vals$DD.CNV.all.obs_p,vals$SCZ.CNV.all.obs_p,
                           vals$DD_SCZ.CNV.all.obs_p,vals$CNCR.CNV.all.obs_p,
                           na.rm=T)
  if(!(all(is.na(c(vals$DD.CNV.all.perm_p,vals$SCZ.CNV.all.perm_p,
                   vals$DD_SCZ.CNV.all.perm_p,vals$CNCR.CNV.all.perm_p))))){
    MIN.CNV.all.perm_p <- min(vals$DD.CNV.all.perm_p,vals$SCZ.CNV.all.perm_p,
                              vals$DD_SCZ.CNV.all.perm_p,vals$CNCR.CNV.all.perm_p,
                              na.rm=T)
  }else{
    MIN.CNV.all.perm_p <- NA
  }
  #All CNVs, noncoding only
  MIN.CNV.noncoding.obs_p <- min(vals$DD.CNV.noncoding.obs_p,vals$SCZ.CNV.noncoding.obs_p,
                           vals$DD_SCZ.CNV.noncoding.obs_p,vals$CNCR.CNV.noncoding.obs_p,
                           na.rm=T)
  if(!(all(is.na(c(vals$DD.CNV.noncoding.perm_p,vals$SCZ.CNV.noncoding.perm_p,
                   vals$DD_SCZ.CNV.noncoding.perm_p,vals$CNCR.CNV.noncoding.perm_p))))){
    MIN.CNV.noncoding.perm_p <- min(vals$DD.CNV.noncoding.perm_p,vals$SCZ.CNV.noncoding.perm_p,
                              vals$DD_SCZ.CNV.noncoding.perm_p,vals$CNCR.CNV.noncoding.perm_p,
                              na.rm=T)
  }else{
    MIN.CNV.noncoding.perm_p <- NA
  }
  #DELs, coding + noncoding
  MIN.DEL.all.obs_p <- min(vals$DD.DEL.all.obs_p,vals$SCZ.DEL.all.obs_p,
                           vals$DD_SCZ.DEL.all.obs_p,vals$CNCR.DEL.all.obs_p,
                           na.rm=T)
  if(!(all(is.na(c(vals$DD.DEL.all.perm_p,vals$SCZ.DEL.all.perm_p,
                   vals$DD_SCZ.DEL.all.perm_p,vals$CNCR.DEL.all.perm_p))))){
    MIN.DEL.all.perm_p <- min(vals$DD.DEL.all.perm_p,vals$SCZ.DEL.all.perm_p,
                              vals$DD_SCZ.DEL.all.perm_p,vals$CNCR.DEL.all.perm_p,
                              na.rm=T)
  }else{
    MIN.DEL.all.perm_p <- NA
  }
  #DELs, noncoding only
  MIN.DEL.noncoding.obs_p <- min(vals$DD.DEL.noncoding.obs_p,vals$SCZ.DEL.noncoding.obs_p,
                                 vals$DD_SCZ.DEL.noncoding.obs_p,vals$CNCR.DEL.noncoding.obs_p,
                                 na.rm=T)
  if(!(all(is.na(c(vals$DD.DEL.noncoding.perm_p,vals$SCZ.DEL.noncoding.perm_p,
                   vals$DD_SCZ.DEL.noncoding.perm_p,vals$CNCR.DEL.noncoding.perm_p))))){
    MIN.DEL.noncoding.perm_p <- min(vals$DD.DEL.noncoding.perm_p,vals$SCZ.DEL.noncoding.perm_p,
                                    vals$DD_SCZ.DEL.noncoding.perm_p,vals$CNCR.DEL.noncoding.perm_p,
                                    na.rm=T)
  }else{
    MIN.DEL.noncoding.perm_p <- NA
  }
  #DUPs, coding + noncoding
  MIN.DUP.all.obs_p <- min(vals$DD.DUP.all.obs_p,vals$SCZ.DUP.all.obs_p,
                           vals$DD_SCZ.DUP.all.obs_p,vals$CNCR.DUP.all.obs_p,
                           na.rm=T)
  if(!(all(is.na(c(vals$DD.DUP.all.perm_p,vals$SCZ.DUP.all.perm_p,
                   vals$DD_SCZ.DUP.all.perm_p,vals$CNCR.DUP.all.perm_p))))){
    MIN.DUP.all.perm_p <- min(vals$DD.DUP.all.perm_p,vals$SCZ.DUP.all.perm_p,
                              vals$DD_SCZ.DUP.all.perm_p,vals$CNCR.DUP.all.perm_p,
                              na.rm=T)
  }else{
    MIN.DUP.all.perm_p <- NA
  }
  #DUPs, noncoding only
  MIN.DUP.noncoding.obs_p <- min(vals$DD.DUP.noncoding.obs_p,vals$SCZ.DUP.noncoding.obs_p,
                                 vals$DD_SCZ.DUP.noncoding.obs_p,vals$CNCR.DUP.noncoding.obs_p,
                                 na.rm=T)
  if(!(all(is.na(c(vals$DD.DUP.noncoding.perm_p,vals$SCZ.DUP.noncoding.perm_p,
                   vals$DD_SCZ.DUP.noncoding.perm_p,vals$CNCR.DUP.noncoding.perm_p))))){
    MIN.DUP.noncoding.perm_p <- min(vals$DD.DUP.noncoding.perm_p,vals$SCZ.DUP.noncoding.perm_p,
                                    vals$DD_SCZ.DUP.noncoding.perm_p,vals$CNCR.DUP.noncoding.perm_p,
                                    na.rm=T)
  }else{
    MIN.DUP.noncoding.perm_p <- NA
  }
  #Summary minimums across classes
  MIN.CNV.obs_p <- min(MIN.CNV.all.obs_p,MIN.CNV.noncoding.obs_p)
  if(!(all(is.na(c(MIN.CNV.all.perm_p,MIN.CNV.noncoding.perm_p))))){
    MIN.CNV.perm_p <- min(MIN.CNV.all.perm_p,MIN.CNV.noncoding.perm_p,na.rm=T)
  }else{
    MIN.CNV.perm_p <- NA
  }
  MIN.DEL.obs_p <- min(MIN.DEL.all.obs_p,MIN.DEL.noncoding.obs_p)
  if(!(all(is.na(c(MIN.DEL.all.perm_p,MIN.DEL.noncoding.perm_p))))){
    MIN.DEL.perm_p <- min(MIN.DEL.all.perm_p,MIN.DEL.noncoding.perm_p,na.rm=T)
  }else{
    MIN.DEL.perm_p <- NA
  }
  MIN.DUP.obs_p <- min(MIN.DUP.all.obs_p,MIN.DUP.noncoding.obs_p)
  if(!(all(is.na(c(MIN.DUP.all.perm_p,MIN.DUP.noncoding.perm_p))))){
    MIN.DUP.perm_p <- min(MIN.DUP.all.perm_p,MIN.DUP.noncoding.perm_p,na.rm=T)
  }else{
    MIN.DUP.perm_p <- NA
  }
  MIN.ANY.obs_p <- min(MIN.CNV.obs_p,MIN.DEL.obs_p,MIN.DUP.obs_p)
  if(!(all(is.na(c(MIN.CNV.perm_p,MIN.DEL.perm_p,MIN.DUP.perm_p))))){
    MIN.ANY.perm_p <- min(MIN.CNV.perm_p,MIN.DEL.perm_p,MIN.DUP.perm_p,na.rm=T)
  }else{
    MIN.ANY.perm_p <- NA
  }
  
  #Return
  return(c(MIN.CNV.all.obs_p,MIN.CNV.all.perm_p,
    MIN.CNV.noncoding.obs_p,MIN.CNV.noncoding.perm_p,
    MIN.DEL.all.obs_p,MIN.DEL.all.perm_p,
    MIN.DEL.noncoding.obs_p,MIN.DEL.noncoding.perm_p,
    MIN.DUP.all.obs_p,MIN.DUP.all.perm_p,
    MIN.DUP.noncoding.obs_p,MIN.DUP.noncoding.perm_p,
    MIN.CNV.obs_p,MIN.CNV.perm_p,
    MIN.DEL.obs_p,MIN.DEL.perm_p,
    MIN.DUP.obs_p,MIN.DUP.perm_p,
    MIN.ANY.obs_p,MIN.ANY.perm_p))
})))
names(merged_results) <- c("MIN.CNV.all.obs_p","MIN.CNV.all.perm_p",
                           "MIN.CNV.noncoding.obs_p","MIN.CNV.noncoding.perm_p",
                           "MIN.DEL.all.obs_p","MIN.DEL.all.perm_p",
                           "MIN.DEL.noncoding.obs_p","MIN.DEL.noncoding.perm_p",
                           "MIN.DUP.all.obs_p","MIN.DUP.all.perm_p",
                           "MIN.DUP.noncoding.obs_p","MIN.DUP.noncoding.perm_p",
                           "MIN.CNV.obs_p","MIN.CNV.perm_p",
                           "MIN.DEL.obs_p","MIN.DEL.perm_p",
                           "MIN.DUP.obs_p","MIN.DUP.perm_p",
                           "MIN.ANY.obs_p","MIN.ANY.perm_p")
results_out <- cbind(x,merged_results)
names(results_out)[1] <- "#chr"
write.table(results_out,args[2],col.names=T,row.names=F,quote=F,sep="\t")


