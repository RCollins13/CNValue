#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate swarm plot of exonic vs whole-gene DEL vs DUP

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

#####Read deletion data
DEL.ex <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DEL.ex.means <- apply(DEL.ex[,-1],1,mean,na.rm=T)
DEL.ex.means.log <- log2(DEL.ex.means)
DEL.wg <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DEL.wg.means <- apply(DEL.wg[,-1],1,mean,na.rm=T)
DEL.wg.means.log <- log2(DEL.wg.means)

#####Read duplication data
DUP.ex <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DUP_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DUP.ex.means <- apply(DUP.ex[,-1],1,mean,na.rm=T)
DUP.ex.means.log <- log2(DUP.ex.means)
DUP.wg <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DUP_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DUP.wg.means <- apply(DUP.wg[,-1],1,mean,na.rm=T)
DUP.wg.means.log <- log2(DUP.wg.means)

#####Prepare plotting area
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/ex_vs_wg_DEL_vs_DUP.swarms.png",sep=""),
    width=4*0.8,height=3,units="in",res=1000)
par(bty="n",mar=c(0.5,2.5,0.5,0.5))
plot(x=c(0.5,4.5),y=log2(c(1/3,7)),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")

#####Draw gridlines
abline(h=log2(c(1/6,1/4,1/2,1,2,4,6)),col=cols.CTRL[3])
abline(h=0,lwd=1.5,col="gray20")

#####Draw y-axis
axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(2,at=log2(c(1/6,1/4,1/2,1,2,4,6)),labels=NA)
axis(2,at=log2(c(1/6,1/4,1/2,1,2,4,6)),tick=F,
     labels=c("1/6","1/4","1/2",1,2,4,6),las=2)

#####Iterate over groups & plot swarms
means <- list(DEL.ex.means.log,DUP.ex.means.log,
              DEL.wg.means.log,DUP.wg.means.log)
mean.cols <- rep(c("red","blue"),2)
sapply(1:4,function(i){
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
  segments(x0=i,x1=i,y0=lower,y1=upper)
  segments(x0=i-0.05,x1=i+0.05,y0=c(lower,upper),y1=c(lower,upper))
  points(x=i,y=median(means[[i]],na.rm=T),
         pch=23,bg="white")
})

#####Draw closing lines
abline(h=par("usr")[3])
# axis(1,at=1:2,tck=F,labels=NA)
# axis(1,at=3:4,tck=F,labels=NA)

#####Close device
dev.off()

#####Mann-Whitney tests
options(scipen=-1000)
wilcox.test(DEL.ex.means.log,
            DUP.ex.means.log)$p.value
wilcox.test(DEL.ex.means.log,
            DEL.wg.means.log)$p.value
wilcox.test(DUP.wg.means.log,
            DEL.wg.means.log,
            alternative="greater")$p.value
wilcox.test(DUP.wg.means.log,
            DUP.ex.means.log,
            alternative="greater")$p.value
options(scipen=1000)







