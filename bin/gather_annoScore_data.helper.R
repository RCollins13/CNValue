#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Script to calculate annoScore data from CNV-anno pairs

##Command line arguments:
# CASE: path to CNV-anno pairs for cases
# CTRL: path to CNV-anno pairs for controls
# IDs: path to list of all element IDs
# OUTDIR: path to output directory

##Local dev testing parameters
CASE <- "~/scratch/tmp.qmLhWlBisx/tmp.BP1Oqq7UV1"
CTRL <- "~/scratch/tmp.qmLhWlBisx/tmp.IRjLHfdNbB"
IDs <- "~/scratch/tmp.qmLhWlBisx/tmp.D8Fqduq2r9"
OUTDIR <- "~/scratch/tmp.qmLhWlBisx/"

#####Set params
options(scipen=1000,stringsAsFactors=F)

#####Load requirements
require("optparse")

#####################################################
####Helper function to calculate all relevant metrics
#####################################################
calcCNVdata <- function(path,IDs){
  #Read CNV-anno pairs
  df <- read.table(path,header=F)

  #Count number of annos per CNV
  annosPerCNV <- as.data.frame(table(df[,1]))

  #Count number of CNVs per anno
  CNVsPerAnno <- as.data.frame(table(df[,2]))
  CNVsPerAnno <- t(sapply(IDs,function(ID){
    if(ID %in% as.character(CNVsPerAnno[,1])){
      counts <- CNVsPerAnno[which(CNVsPerAnno[,1]==ID),2]
    }else{
      counts <- 0
    }
    return(c(ID,counts))
  }))

  #Join CNV-anno pairs with number of annos per CNV
  merged.df <- merge(df,annosPerCNV,by=1)
  merged.df[,3] <- 1/merged.df[,3]

  #Get sum of weighted CNVs per anno
  weightedCNVsPerAnno <- as.data.frame(t(sapply(IDs,function(ID){
    return(c(as.character(ID),
             sum(merged.df[which(merged.df[,2]==as.character(ID)),3])))
  })))

  #Return list of results
  return(list(annosPerCNV,CNVsPerAnno,weightedCNVsPerAnno))
}

#######################################
#####Rscript command line functionality
#######################################
#List of Rscript options
option_list <- list()

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] CASE CTRL OUTDIR",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#Assigns positional arguments to named variables
CASE <- as.character(args$args[1])
CTRL <- as.character(args$args[2])
OUTDIR <- as.character(args$args[3])

#Checks for appropriate positional arguments
if(length(args$args) != 4) {
  print_help(parser)
  stop("Must provide four positional arguments")
}

#Read list of IDs and convert to vector
IDs <- as.vector(as.character(read.table(IDs,header=F)[,1]))

#Process CASE data
CASE.dat <- calcCNVdata(CASE,IDs)

#Process CTRL data
CTRL.dat <- calcCNVdata(CTRL,IDs)

#Write CASE outputs
write.table(CASE.dat[[1]],
            paste(OUTDIR,"/CASE.annosPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CASE.dat[[2]],
            paste(OUTDIR,"/CASE.CNVsPerAnno.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CASE.dat[[3]],
            paste(OUTDIR,"/CASE.weightedCNVsPerAnno.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)

#Write CTRL outputs
write.table(CTRL.dat[[1]],
            paste(OUTDIR,"/CTRL.annosPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CTRL.dat[[2]],
            paste(OUTDIR,"/CTRL.CNVsPerAnno.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CTRL.dat[[3]],
            paste(OUTDIR,"/CTRL.weightedCNVsPerAnno.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)

