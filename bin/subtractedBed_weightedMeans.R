#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Calculates mean over a bed interval (a la bedtools coverage) but allows for
#intervals to be subtracted from original interval

#Input:
# subtracted bed file
# output file

#Function to calcuate weighted means
wmeans <- function(input,output){
	options(stringsAsFactors=F)
	x <- read.table(input,header=F)
	weighted <- sapply(unique(x[,4]),function(bin){
	  splits <- x[which(x[,4]==bin),]
	  sizes <- splits[,3]-splits[,2]
	  return(sum(sizes*splits[,5])/sum(sizes))
	})
	weighted <- data.frame(names(weighted),weighted)
	allbins <- read.table("/data/talkowski/rlc47/TAD_intolerance/data/annotations/GRCh37.autosomes.master_bins.bed.gz",header=F)
	res <- merge(allbins,weighted,sort=F,by.x=4,by.y=1,all.x=T)[,c(2:4,1,5)]
	write.table(res,output,col.names=F,row.names=F,sep="\t",quote=F)
}

#Read args
args <- commandArgs(trailingOnly=T)

#Run function
wmeans(args[1],args[2])