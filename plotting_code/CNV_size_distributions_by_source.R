#!/usr/bin/env R

#Copyright (c) 2017 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#Plot size distributions of CNVs by source

####Set parameters####
options(scipen=6)
# sources <- c("SSC_CTRL","Shaikh_CTRL","Suktitipat_CTRL","Cooper_CTRL","Uddin_CTRL","Vogler_CTRL","Coe_CTRL",
#              "PGC_CTRL","TCGA_CTRL","TCGA_CNCR","PGC_GERM","Coe_GERM","Talkowski_GERM","SSC_GERM")
sources <- c("SSC_CTRL","Shaikh_CTRL","Suktitipat_CTRL","Cooper_CTRL","Uddin_CTRL","Vogler_CTRL",
             "PGC_CTRL","PGC_GERM","Coe_GERM","Talkowski_GERM","SSC_GERM","TCGA_CNCR")
source.col <- c(rep("#767071",7),rep("#04B150",4),"#F5717D",
                "#767071","#04B150","#F5717D","white")
WRKDIR <- "/Users/collins/Desktop/RCollins/Talkowski_Local/CNV_DB/rCNV_map/"
PLOTDIR <- paste(WRKDIR,"plots/",sep="")
#Create PLOTDIR if it doesn't already exist
if(!dir.exists(PLOTDIR)){
  dir.create(PLOTDIR)
}

####Load libraries####
require(vioplot)
library(plotrix)

####Function to plot all violins per rCNV type####
plotCNV <- function(CNV="CNV"){
  pdf(paste(PLOTDIR,CNV,"_size_by_source.pdf",sep=""),height=4,width=8)
  #Read data
  CNV.sizes <- sapply(sources,function(source){
    data <- as.vector(read.table(paste(WRKDIR,"plot_data/sizes_by_source/",CNV,"/",source,"_",CNV,".sizes.txt",sep=""),header=F)[,1])
  })
  #Add master CTRL, GERM, CNCR, and ALL vectors
  CNV.sizes$CTRL <- c(unlist(sapply(1:7,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$GERM <- c(unlist(sapply(8:11,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$CNCR <- c(unlist(sapply(12,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$ALL <- c(unlist(sapply(1:12,function(i){return(CNV.sizes[[i]])})))
  #Set log scaling
  logvect <- log10(as.vector((sapply(4:7,function(x){return(1:9*(10^x))}))))
  #Prep plot
  par(mar=c(4.1,4.1,2.6,0.6))
  plot(x=c(0,length(CNV.sizes)),y=log10(c(50000,5000000)),
       type="n",ylab="",yaxt="n",xlab="",xaxt="n",yaxs="i",xaxs="i",
       main=paste(CNV,"Size (log10-scaled)",sep=" "))
  abline(h=logvect,col="gray90",lwd=0.5)
  abline(h=logvect[seq(5,39,9)],col="gray80")
  abline(h=logvect[seq(1,39,9)],col="gray70",lwd=1.5)
  axis(2,at=log10(c(20000,50000,100000,500000,1000000,5000000)),
       labels=c("20kb","50kb","100kb","500kb","1Mb","5Mb"),las=2)
  #Plot violins
  sapply(1:length(CNV.sizes),function(i){
    vioplot(log10(CNV.sizes[[i]]),h=1/75,names=NA,col=source.col[i],at=i-0.5,add=T,pchMed=18)
  })
  #Add xaxis
  axis(1,at=0.5:15.5,tick=F,line=-0.8,labels=c(sources,"CTRL_ALL","GERM_ALL","CNCR_ALL","ALL"),
       las=2,cex.axis=0.7)
  abline(v=12,lwd=3)
  dev.off()
}

#####Run####
plotCNV("CNV")
plotCNV("DEL")
plotCNV("DUP")


####Function to plot all violins per urCNV type####
ploturCNV <- function(CNV="CNV"){
  pdf(paste(PLOTDIR,CNV,"_urCNV_size_by_source.pdf",sep=""),height=4,width=8)
  #Read data
  CNV.sizes <- sapply(sources,function(source){
    data <- as.vector(read.table(paste(WRKDIR,"plot_data/sizes_by_source/",CNV,"_urCNV/",source,"_",CNV,".sizes.txt",sep=""),header=F)[,1])
  })
  #Add master CTRL, GERM, CNCR, and ALL vectors
  CNV.sizes$CTRL <- c(unlist(sapply(1:7,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$GERM <- c(unlist(sapply(8:11,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$CNCR <- c(unlist(sapply(12,function(i){return(CNV.sizes[[i]])})))
  CNV.sizes$ALL <- c(unlist(sapply(1:12,function(i){return(CNV.sizes[[i]])})))
  #Set log scaling
  logvect <- log10(as.vector((sapply(4:7,function(x){return(1:9*(10^x))}))))
  #Prep plot
  par(mar=c(4.1,4.1,2.6,0.6))
  plot(x=c(0,length(CNV.sizes)),y=log10(c(50000,5000000)),
       type="n",ylab="",yaxt="n",xlab="",xaxt="n",yaxs="i",xaxs="i",
       main=paste("ur",CNV," Size (log10-scaled)",sep=""))
  abline(h=logvect,col="gray90",lwd=0.5)
  abline(h=logvect[seq(5,39,9)],col="gray80")
  abline(h=logvect[seq(1,39,9)],col="gray70",lwd=1.5)
  axis(2,at=log10(c(20000,50000,100000,500000,1000000,5000000)),
       labels=c("20kb","50kb","100kb","500kb","1Mb","5Mb"),las=2)
  #Plot violins
  sapply(1:length(CNV.sizes),function(i){
    vioplot(log10(CNV.sizes[[i]]),h=1/75,names=NA,col=source.col[i],at=i-0.5,add=T,pchMed=18)
  })
  #Add xaxis
  axis(1,at=0.5:15.5,tick=F,line=-0.8,labels=c(sources,"CTRL_ALL","GERM_ALL","CNCR_ALL","ALL"),
       las=2,cex.axis=0.7)
  abline(v=12,lwd=3)
  dev.off()
}

#####Run####
ploturCNV("CNV")
ploturCNV("DEL")
ploturCNV("DUP")
