#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate swarm plot of genic burden ORs by phenotype for Fig3a

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(beeswarm)
require(plotrix)

#####Read data
OR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                       "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
#Average across all NEURO, SOMA, and CNCR
GERM.mean <- log2(apply(OR[,2:23],1,mean,na.rm=T))
NEURO.mean <- log2(apply(OR[,2:13],1,mean,na.rm=T))
SOMA.mean <- log2(apply(OR[,14:23],1,mean,na.rm=T))
CNCR.mean <- log2(apply(OR[,24:36],1,mean,na.rm=T))
ALL.mean <- log2(apply(OR[,-1],1,mean,na.rm=T))

#####Simulate expected gene set if normally distributed around OR=1
EXP.mean <- rnorm(mean=0,sd=sd(ALL.mean,na.rm=T),n=length(ALL.mean))

#####Prepare plotting area
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/genic_ORs_byPheno.swarms.png",sep=""),
    width=4,height=3,units="in",res=1000)
par(bty="n",mar=c(0.5,2.5,0.5,0.5))
plot(x=c(0.5,5.5),y=log2(c(1/3,6)),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")

#####Draw gridlines
abline(h=log2(c(1/6,1/4,1/2,1,2,4,6)),col=cols.CTRL[4])
abline(h=0,lwd=1.5,col="gray20")

#####Draw y-axis
axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(2,at=log2(c(1/6,1/4,1/2,1,2,4,6)),labels=NA)
axis(2,at=log2(c(1/6,1/4,1/2,1,2,4,6)),tick=F,
     labels=c("1/6","1/4","1/2",1,2,4,6),las=2)

#####Iterate over groups & plot swarms
means <- list(EXP.mean,GERM.mean,NEURO.mean,SOMA.mean,CNCR.mean)
mean.cols <- c(cols.CTRL[1],cols.GERM[1],cols.NEURO[1],cols.SOMA[1],cols.CNCR[1])
sapply(1:5,function(i){
  #Swarm
  beeswarm(means[[i]],add=T,at=i,method="swarm",corral="random",
           pch=19,col=mean.cols[i],cex=0.4)
  #Median & 95% CI
  lower <- median(means[[i]],na.rm=T)-1.96*std.error(means[[i]],na.rm=T)
  upper <- median(means[[i]],na.rm=T)+1.96*std.error(means[[i]],na.rm=T)
  # segments(x0=i-0.1,x1=i+0.1,
  #          y0=median(means[[i]],na.rm=T),
  #          y1=median(means[[i]],na.rm=T),
  #          lwd=1.5)
  points(x=i,y=median(means[[i]],na.rm=T),
         pch=15,cex=0.75)
  segments(x0=i,x1=i,y0=lower,y1=upper)
  segments(x0=i-0.05,x1=i+0.05,y0=c(lower,upper),y1=c(lower,upper))
})

#####Draw closing lines
abline(h=par("usr")[3:4])

#####Close device
dev.off()








