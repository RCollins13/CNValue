#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate histogram of phenotype assignments per subject

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)

#####Read data
df <- read.table(paste(WRKDIR,"plot_data/rCNVdbSummaryFigure/phenotypes_per_sample_hist.txt",sep=""))
df[9,2] <- sum(df[9:nrow(df),2])
df <- df[1:9,]

#####Prepare plot
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/phenotypes_per_subject.histogram.pdf",sep=""),
    width=3,height=2.7)
par(mar=c(1.5,3,0.5,0.5),bty="n")
plot(x=c(-0.1,nrow(df)),y=c(0,1.02*max(df[,2])),
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

#####Plot barplot
rect(xleft=df[,1]-0.9,xright=df[,1]-0.1,
     ybottom=rep(0,times=nrow(df)),ytop=df[,2],
     col="#A5A6A7",lwd=1.75)

#####Add axes
axis(1,at=0.5:7.5,tick=F,line=-0.8,labels=1:8)
axis(1,at=8.5,tick=F,line=-0.8,labels=">8")
axis(2,at=seq(0,35000,2500),col="#A5A6A7",tck=-0.01,labels=NA)
axis(2,at=seq(0,35000,5000),col="#373738",tck=-0.015,labels=NA)
# axis(2,at=seq(5000,25000,10000),col.lab="#373738",tick=F,
#      labels=c("5k","15k","25k"),cex.axis=0.8,las=2,line=-0.4)
axis(2,at=seq(0,30000,10000),tck=-0.02,labels=NA)
axis(2,at=seq(0,30000,5000),tick=F,labels=c(0,paste(seq(5,30,5),"k",sep="")),las=2,line=-0.2)

#####Close plot
dev.off()
