#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Script to filter master bins burden file

#Set params
options(scipen=1000,stringsAsFactors=F)
diseases <- c("DD","SCZ","DD_SCZ","CNCR","ANY_DISEASE")
CNVs <- c("CNV","DEL","DUP","ANY_CNV")
filts <- c("all","coding","noncoding","ANY_FILTER")
actions <- c("include","exclude")

#Load requirements
require("optparse")

#Function to filter
filter.table <- function(df,
                         terms,
                         threshold=0.05,
                         action="include"){
  #Sanity check action
  if(!(action %in% actions)){
    stop("Argument 'action' must be either 'include' or 'exclude'")
  }

  #Get corresponding column indexes
  cidx <- which(names(df) %in% terms)
  if(length(cidx)==0){
    stop("No column headers match your specified terms. Please check your input terms list.")
  }

  #Get rows passing filter for each column
  if(action=="include"){
    pass <- lapply(cidx,function(i){
      return(which(df[,i] <= threshold & !(is.na(df[,i]))))
    })
  }else{
    pass <- lapply(cidx,function(i){
      return(which(df[,i] > threshold & !(is.na(df[,i]))))
    })
  }

  #Find intersection of all passing rows
  union <- pass[[1]]
  if(length(pass)>1){
    for(i in 2:length(pass)){
      union <- intersect(union,pass[[i]])
    }
  }

  #Filter df
  results <- df[union,]

  #Return results
  names(results)[1] <- "#chr"
  return(results)
}

#list of Rscript options
option_list <- list(
  make_option(c("-a", "--action"), type="character", default="include",
              help="include or exclude bins based on criteria [default '%default']",
              metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.05,
              help="p-value threshold for including or excluding bins [default '%default']",
              metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="/dev/stdout",
              help="write output to file [default stdout]",
              metavar="character")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] burden_file terms_list",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#Checks for appropriate positional arguments
if(length(args$args) != 2) {
  print_help(parser)
  stop("Must specify an input burden file and a terms list")
}

#Loads data frame
df <- read.table(args$args[1],comment.char="",header=T)

#Loads terms & converts to vector
terms <- as.vector(read.table(args$args[2],header=F)[,1])


#Filters
output <- filter.table(df,
                       terms,
                       threshold=opts$threshold,
                       action=opts$action)
names(output)[1] <- "#chr"

#Write out
if(opts$outfile != "/dev/stdout"){
  write.table(output,opts$outfile,col.names=T,row.names=F,sep="\t",quote=F)
}else{
  print(output,row.names=F)
}

