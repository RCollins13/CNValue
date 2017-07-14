#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate CNV Z-score vs ExAC LOF Z-score scatterplots

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

#####Helper function to read data, & match with ExAC Z-scores (based on gene symbol)
#Read ExAC Z-scores
ExAC.Z <- read.table(paste(WRKDIR,"plot_data/figure4/ExAC_LoF_Zscores.txt",sep=""),
                     header=T)
#Read CNV data & match w/ExAC Z-scores
readData <- function(pheno,VF,context){
  lapply(list("CNV","DEL","DUP"),function(CNV){
    #Read table
    dat <- read.table(paste(WRKDIR,"plot_data/figure4/geneScore_results/",
                            pheno,"_",CNV,"_",VF,"_",context,
                            ".geneScore_stats.txt",sep=""),
                      comment.char="",header=T)
    #Replace first column header
    names(dat)[1] <- "gene"
    #Match w/ExAC
    dat <- merge(dat,ExAC.Z)
    #Return data
    return(dat)
  })
}




