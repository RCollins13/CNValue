#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate dummy plots for rCNVmap per-gene figure panel A (geneScore schematic)

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#Load libraries
require(MASS)

#################################################
#####Exponential decay model (for weighting CNVs)
#################################################
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/CNV_weighting_expDecay.pdf",sep=""),
    width=2,height=2)
par(mar=c(0.3,0.3,0.2,0.2))
plot(x=c(0.93,6.07),y=c(0,1.02),type="n",
     xaxs="i",xaxt="n",xlab="",
     yaxs="i",yaxt="n",ylab="")
abline(h=seq(0,1,0.25),v=1:10,col=cols.CTRL[3])
polygon(x=c(seq(0,10,0.1),rev(seq(0,10,0.1))),
        y=c(10/seq(0,100,1),rep(0,101)),
        border=NA,col=cols.CTRL[1])
points(x=seq(0,10,0.1),
       y=10/seq(0,100,1),
       type="l",lwd=3)
points(x=1:6,y=c(1/(1:6)),
       pch=21,lwd=3,bg="white")
axis(1,at=1:6,labels=NA)
axis(2,at=seq(0,1,0.2),labels=NA)
dev.off()

########################################
#####Regression model (for fitting CNVs)
########################################
par(mar=rep(0.2,4))
plot(x=c(-3,3),y=c(-3,3),type="n",
     xaxs="i",xaxt="n",xlab="",
     yaxs="i",yaxt="n",ylab="")
points(rnorm(1000),
       rnorm(1000))









