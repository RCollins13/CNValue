#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to run statistical comparisons of full CNV dataset for paper & figure legends

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)

#####Mann-Whitney U for CNV size differences
#germline cases vs controls, all CNVs
cases <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.GERM.E2.all.txt",sep=""))[,1]
controls <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.CTRL.E2.all.txt",sep=""))[,1]
mean(cases)/mean(controls)
wilcox.test(cases,controls,alternative="greater")$p.value
#cases vs controls, coding CNVs
cases <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.GERM.E2.coding.txt",sep=""))[,1]
cancer <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.CNCR.E2.coding.txt",sep=""))[,1]
controls <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.CTRL.E2.coding.txt",sep=""))[,1]
mean(cases)/mean(controls)
format(wilcox.test(cases,controls,alternative="greater")$p.value,scientific=T)
mean(cancer)/mean(controls)
format(wilcox.test(cancer,controls,alternative="greater")$p.value,scientific=T)
#cases vs controls, noncoding CNVs
cases <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.GERM.E2.noncoding.txt",sep=""))[,1]
cancer <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.CNCR.E2.noncoding.txt",sep=""))[,1]
controls <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.CTRL.E2.noncoding.txt",sep=""))[,1]
mean(cases)/mean(controls)
format(wilcox.test(cases,controls,alternative="greater")$p.value,scientific=T)
mean(cancer)/mean(controls)
format(wilcox.test(cancer,controls,alternative="greater")$p.value,scientific=T)


