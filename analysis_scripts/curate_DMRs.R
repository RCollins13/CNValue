#!/usr/bin/env Rscript

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Split DMRs by donor source from Ziller et al., Nature, 2013

#Master function
splitDMRs <- function(path,outdir){
  #Read data & subset to coordinates + sample means 
  dat <- read.table(path,header=T)[,c(2:4,6:29)]
  
  #Iterate over all DMRs and return the column of the 
  #most hypomethylated sample
  mins <- as.vector(unlist(apply(dat[,-c(1:3)],1,function(vals){
    return(which(vals==min(vals,na.rm=T))[1]+3)
  })))
  
  #Iterate over all samples and write tissue-specific 
  #hypomethylated sites to one bed file per sample
  sapply(4:ncol(dat),function(i){
    DMRs <- dat[which(mins==i),1:3]
    DMRs[,1] <- sub("chr","",DMRs[,1])
    write.table(DMRs,paste(outdir,"/",names(dat)[i],".DMRs.bed",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
  })
}

#Read path & outdir from command line, run split
args <- commandArgs(trailingOnly=T)
splitDMRs(args[1],args[2])