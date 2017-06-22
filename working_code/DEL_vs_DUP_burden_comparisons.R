#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Exploratory code for del/dup contrasts from gene & anno burden tests

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Function to load & clean data
load.results <- function(CNV,VF,filt){
  #Load OR point estimates
  OR <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                         CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
  #log2-transform ORs & remove Inf
  OR[,-1] <- log2(OR[,-1])
  is.na(OR[,-1]) <- sapply(OR[,-1],is.infinite)
  #Calculate average OR for GERM, NEURO, SOMA, CNCR
  OR$GERM.mean <- apply(OR[,4:23],1,mean,na.rm=T)
  OR$NEURO.mean <- apply(OR[,4:13],1,mean,na.rm=T)
  OR$SOMA.mean <- apply(OR[,14:23],1,mean,na.rm=T)
  OR$CNCR.mean <- apply(OR[,24:36],1,mean,na.rm=T)

  #Load p-vals
  p <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                         CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)
  #-log10 transform p values
  p[,-1] <- apply(p[,-1],2,as.numeric)
  p[,-1] <- -log10(p[,-1])
  p[,-1] <- apply(p[,-1],2,function(vals){
    vals[which(is.infinite(vals))] <- 300
    return(vals)
  })
  #Calculate average p for GERM, NEURO, SOMA, CNCR
  p$GERM.mean <- apply(p[,4:23],1,mean,na.rm=T)
  p$NEURO.mean <- apply(p[,4:13],1,mean,na.rm=T)
  p$SOMA.mean <- apply(p[,14:23],1,mean,na.rm=T)
  p$CNCR.mean <- apply(p[,24:36],1,mean,na.rm=T)

  #Return data
  return(list(OR[,c(1,37:40)],
              p[,c(1,37:40)]))
}

#####Load deletion data
DEL <- load.results("DEL","E2","noncoding")
DUP <- load.results("DUP","E2","noncoding")

#####Test scatterplot
par(mfrow=c(2,2))
sapply(c("E2","E3","E4","N1"),function(VF){
  DEL <- load.results("DEL",VF,"noncoding")
  DUP <- load.results("DUP",VF,"noncoding")
  plot(DEL[[1]]$CNCR.mean,
       DUP[[1]]$CNCR.mean,
       xlim=c(-4,4),ylim=c(-4,4),
       pch=19,col=adjustcolor("black",alpha=0.3))
  abline(h=0,v=0)
  abline(lm(DUP[[1]]$CNCR.mean ~ DEL[[1]]$CNCR.mean))
  cor(DEL[[1]]$CNCR.mean,
      DUP[[1]]$CNCR.mean,
      use="complete.obs")
  cor.test(DEL[[1]]$CNCR.mean,
            DUP[[1]]$CNCR.mean)$p.value
})

# #Assign colors based on significance
# colvect <- sapply(1:length(DEL[[1]]$GERM.mean),function(i){
#   DEL.p <- DEL[[2]]$GERM.mean[i]
#   if(is.nan(DEL.p) | is.na(DEL.p)){
#     DEL.p <- 0
#   }
#   DUP.p <- DUP[[2]]$GERM.mean[i]
#   if(is.nan(DUP.p) | is.na(DUP.p)){
#     DUP.p <- 0
#   }
#   min.p <- 10^-(0.05/length(DEL[[1]]$GERM.mean))
#   if(DEL.p>=min.p){
#     if(DUP.p>=min.p){
#       return(cols.GERM[1])
#     }else{
#       return("red")
#     }
#   }else{
#     if(DUP.p>=min.p){
#       return("blue")
#     }else{
#       return("black")
#     }
#   }
# })











