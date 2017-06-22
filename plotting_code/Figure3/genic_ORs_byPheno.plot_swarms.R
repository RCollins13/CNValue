#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate swarm plots and swarm plot matrices of annotation set CNV burden results

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

#####Load libraries
require(beeswarm)

#####Read data
OR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                       "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
#Average across all NEURO, SOMA, and CNCR
GERM.mean <- apply(OR[,4:23],1,mean,na.rm=T)
NEURO.mean <- apply(OR[,4:13],1,mean,na.rm=T)
SOMA.mean <- apply(OR[,14:23],1,mean,na.rm=T)
CNCR.mean <- apply(OR[,24:36],1,mean,na.rm=T)

#####Prepare plotting area
plot(x=c(0,6),y=log2(c(1/6,6)),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")










