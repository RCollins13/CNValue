#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Calculate differences between CNV sizes per cohort

#Set params
options(scipen=1000,stringsAsFactors=F)
WRKDIR <- "/Users/collins/Desktop/RCollins/Talkowski_Local/CNV_DB/rCNV_map/"

#Read data
CTRL.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/CTRL.DEL.sizes.txt",sep=""),header=F)[,1])
CTRL.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/CTRL.DUP.sizes.txt",sep=""),header=F)[,1])
CTRL.CNV <- c(CTRL.DEL,CTRL.DUP)
DD.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/DD.DEL.sizes.txt",sep=""),header=F)[,1])
DD.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/DD.DUP.sizes.txt",sep=""),header=F)[,1])
DD.CNV <- c(DD.DEL,DD.DUP)
SCZ.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/SCZ.DEL.sizes.txt",sep=""),header=F)[,1])
SCZ.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/SCZ.DUP.sizes.txt",sep=""),header=F)[,1])
SCZ.CNV <- c(SCZ.DEL,SCZ.DUP)
CNCR.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/CNCR.DEL.sizes.txt",sep=""),header=F)[,1])
CNCR.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/CNV_sizes/CNCR.DUP.sizes.txt",sep=""),header=F)[,1])
CNCR.CNV <- c(CNCR.DEL,CNCR.DUP)
AFFECTED.DEL <- c(DD.DEL,SCZ.DEL,CNCR.DEL)
AFFECTED.DUP <- c(DD.DUP,SCZ.DUP,CNCR.DUP)
AFFECTED.CNV <- c(AFFECTED.DEL,AFFECTED.DUP)
GERM.DEL <- c(DD.DEL,SCZ.DEL,CTRL.DEL)
GERM.DUP <- c(DD.DUP,SCZ.DUP,CTRL.DUP)
GERM.CNV <- c(GERM.DEL,GERM.DUP)
GERM_AFF.DEL <- c(DD.DEL,SCZ.DEL)
GERM_AFF.DUP <- c(DD.DUP,SCZ.DUP)
GERM_AFF.CNV <- c(GERM_AFF.DEL,GERM_AFF.DUP)

#CNCR vs GERM
t.test(CNCR.CNV,GERM.CNV,alternative="greater")$p.value
mean(CNCR.CNV)/mean(GERM.CNV)

#DD vs CTRL
t.test(DD.CNV,CTRL.CNV,alternative="greater")$p.value
mean(DD.CNV)/mean(CTRL.CNV)

#SCZ vs CTRL - CNV
t.test(SCZ.CNV,CTRL.CNV,alternative="greater")$p.value
mean(SCZ.CNV)/mean(CTRL.CNV)

#SCZ vs CTRL - DEL
t.test(SCZ.DEL,CTRL.DEL,alternative="greater")$p.value
mean(SCZ.DEL)/mean(CTRL.DEL)

#SCZ vs CTRL - DUP
t.test(SCZ.DUP,CTRL.DUP,alternative="greater")$p.value
mean(SCZ.DUP)/mean(CTRL.DUP)








