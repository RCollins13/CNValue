#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

# Genome-wide Fisher's exact test of CNV burdens in cases vs controls

#Helper function to read & format data frame
read.CNVbed <- function(path){
  #Check if file exists
  if(!file.exists(path)){
    stop(paste("Input file ",path," not found",sep=""))
  }
  #Read data and cut to first three/four and last columns
  dat <- read.table(path,header=F,sep="\t")
  if(ncol(dat)>4){
    dat <- dat[,c(1:4,ncol(dat))]
    colnames(dat) <- c("Chr","Start","End","ID","CNVs")
  }else{
    colnames(dat) <- c("Chr","Start","End","CNVs")
  }
  #Return cleaned data frame
  return(dat)
}

#Helper function to merge case & control data frames
merge.case.control <- function(case,control){
  #Truncate data frames if inconsistent between cases and controls
  if(ncol(case) != ncol(control)){
    case <- case[,c(1:3,ncol(case))]
    control <- control[,c(1:3,ncol(case))]
  }
  #Merge case & control
  if(ncol(case)>4){
    dat <- merge(case,control,by=1:4,sort=F)
  }else{
    dat <- merge(case,control,by=1:3,sort=F)
  }
  if(nrow(dat)!=nrow(case) | nrow(dat)!=nrow(control)){
    warning("Bins not identical between case and control input files")
  }
  #Rename columns
  if(ncol(dat)==6){
    names(dat) <- c("Chr","Start","End","ID","Case","Control")
  }else{
    names(dat) <- c("Chr","Start","End","Case","Control")
  }
  #Return merged df
  return(dat)
}

#Helper function to run Fisher's Exact Tests on all bins with at least one case CNV
fisher.bins <- function(df){
  #Get column indexes for cases and controls
  case.idx <- ncol(df)-1
  control.idx <- ncol(df)
  #Instantiate p-value and OR vectors
  df$Fisher.p <- NA
  df$OR <- NA
  #Get medians of cases and controls
  case.med <- median(df[,case.idx],na.rm=T)
  control.med <- median(df[,control.idx],na.rm=T)
  #Stop calculations if cases or controls have median of zero CNVs per bin
  if(case.med==0 | control.med==0){
    stop("Cases or controls have median of zero CNVs per bin")
  }
  #Calculate odds ratios per bin if both case & control medians > 0
  if(case.med>0 & control.med>0){
    df$OR <- sapply(1:nrow(df),function(i){
      OR <- (df[i,case.idx]/case.med)/(df[i,control.idx]/control.med)
      return(OR)
    })
  }
  #Get list of bins where case CNV count > 0
  test.bins.idx <- which(df[,case.idx]>0)
  #Run Fisher's exact on all bins with > 0 case CNVs
  pvals.nom <- as.vector(unlist(sapply(test.bins.idx,function(i){
    mat <- data.frame(c(control.med,case.med),
                      c(df[i,control.idx],df[i,case.idx]))
    p <- fisher.test(mat,alternative="greater")$p.value
    return(p)
  })))
  df[test.bins.idx,]$Fisher.p <- pvals.nom
  #Return data frame
  return(df)
}

#Master function
runFisherGenome <- function(case.path,control.path,outfile,gzip=T){
  #Read data
  case <- read.CNVbed(case.path)
  control <- read.CNVbed(control.path)
  #Merge data
  merged.df <- merge.case.control(case,control)
  #Run Fisher
  results <- fisher.bins(merged.df)
  #Write results to file
  write.table(results,outfile,col.names=T,row.names=F,sep="\t",quote=F)
  #Gzip (if optioned)
  if(gzip==T){
    system(paste("gzip -f ",outfile,sep=""),wait=T,intern=F)
  }
}

####################################################
####Rscript functionality for command line usage####
####################################################

#Disables factor default
options(stringsAsFactors=F)

#Requires optparse
require(optparse)

#list of Rscript options
option_list <- list(
  make_option(c("-z", "--gzip"),action="store_true",type="logical",default=FALSE,
              help="gzip output file [default '%default']")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] cases.bed controls.bed outfile",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#checks for appropriate positional arguments
if(length(args$args) != 3) {
  print_help(parser)
  stop("Incorrect number of required positional arguments")
}

#Removes trailing ".gz" from output file if gzip==T
if(opts$gzip==T){
  args$args[3] <- unlist(strsplit(args$args[3],split=".gz",fixed=T))
}

#RunFisherGenome
runFisherGenome(case.path=args$args[1],
                control.path=args$args[2],
                outfile=args$args[3],
                gzip=opts$gzip)
