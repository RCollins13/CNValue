#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate ExAC LoF z-scores for genes in significant gene sets for Fig3

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
require(vioplot)

#####Read data
dat <- lapply(list("allGenes","noPhenos","anyPhenos"),function(category){
  dat <- read.table(paste(WRKDIR,"plot_data/figure3/ExAC_LoF_z.",category,".txt",sep=""),
                         header=F)[,1]
  return(dat)
})

#####Prepare plotting area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/geneSet_ExAC_z_scores.boxplots.pdf",sep=""),
    width=1.5,height=5)
par(bty="n",mar=c(1,2.5,0.5,0.5))
plot(x=c(0.5,3.5),y=c(-3,7),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")

#####Add boxplots
boxplot(dat,add=T,notch=T,outline=F,lty=1,staplewex=0,
        xaxt="n",yaxt="n",col=cols.CTRL[2])

#####Add axes
axis(1,at=1:3,labels=NA)
axis(2,at=-6:8,labels=NA,tck=-0.013,col=cols.CTRL[1])
axis(2,at=seq(-6,8,2),labels=NA)
axis(2,at=seq(-6,8,2),tick=F,las=2,line=-0.3)
mtext(2,text="LoF Constraint Z-Score",line=1.3)

#Close device
dev.off()


