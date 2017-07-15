#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Script to calculate geneScore data from CNV-gene pairs

##Command line arguments:
# CASE: path to CNV-gene pairs for cases
# CTRL: path to CNV-gene pairs for controls
# OUTDIR: path to output directory

##Local dev testing parameters
CASE <- "~/scratch/tmp.AK1hUoqv10/tmp.syaIoQzd5g"
CTRL <- "~/scratch/tmp.AK1hUoqv10/tmp.Ev8ObL7KXr"
OUTDIR <- "~/scratch/tmp.AK1hUoqv10/"
#
# #Test run
# df <- adjustCounts(readGeneScores(infile))
# fisher.results <- calcFisherStats(df,nCTRL,nCASE)
# ratio.results <- calcRatioStats(df,nCTRL,nCASE)
# all.results <- cbind(fisher.results,ratio.results[,-c(1:16)])

#####Set params
options(scipen=1000,stringsAsFactors=F)

#####Load requirements
require("optparse")

#####################################################
####Helper function to calculate all relevant metrics
#####################################################
calcCNVdata <- function(path){
  #Read CNV-gene pairs
  df <- read.table(path,header=F)

  #Count number of genes per CNV
  genesPerCNV <- as.data.frame(table(df[,1]))

  #Count number of CNVs per gene
  CNVsPerGene <- as.data.frame(table(df[,2]))

  #Join CNV-gene pairs with number of genes per CNV
  merged.df <- merge(df,genesPerCNV,by=1)
  merged.df[,3] <- 1/merged.df[,3]

  #Get sum of weighted CNVs per gene
  weightedGenesPerCNV <- as.data.frame(t(sapply(CNVsPerGene[,1],function(gene){
    return(c(as.character(gene),
             sum(merged.df[which(merged.df[,2]==as.character(gene)),3])))
  })))

  #Return list of results
  return(list(genesPerCNV,CNVsPerGene,weightedGenesPerCNV))
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
if(length(args$args) != 3) {
  print_help(parser)
  stop("Must provide three positional arguments")
}

#Process CASE data
CASE.dat <- calcCNVdata(CASE)

#Process CTRL data
CTRL.dat <- calcCNVdata(CTRL)

#Write CASE outputs
write.table(CASE.dat[[1]],
            paste(OUTDIR,"/CASE.genesPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CASE.dat[[2]],
            paste(OUTDIR,"/CASE.CNVsPerGene.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CASE.dat[[3]],
            paste(OUTDIR,"/CASE.weightedGenesPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)

#Write CTRL outputs
write.table(CTRL.dat[[1]],
            paste(OUTDIR,"/CTRL.genesPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CTRL.dat[[2]],
            paste(OUTDIR,"/CTRL.CNVsPerGene.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)
write.table(CTRL.dat[[3]],
            paste(OUTDIR,"/CTRL.weightedGenesPerCNV.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)

