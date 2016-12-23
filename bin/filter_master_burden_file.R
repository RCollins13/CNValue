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
filts <- c("all","noncoding","either_filter")
actions <- c("include","exclude")

#Load requirements
require("optparse")

#Function to filter
filter.table <- function(df,
                         disease="ANY_DISEASE",
                         CNV="ANY_CNV",
                         filt="either_filter",
                         action="include"){
  #Sanity check variables
  if(!(disease %in% diseases)){
    stop("Argument 'disease' must be either 'DD', 'SCZ', 'DD_SCZ', 'CNCR', or 'ANY_DISEASE'")
  }
  if(!(CNV %in% CNVs)){
    stop("Argument 'CNV' must be either 'CNV', 'DEL', 'DUP', or 'ANY_CNV'")
  }
  if(!(filt %in% filts)){
    stop("Argument 'filt' must be either 'all', 'noncoding', or 'either_filter'")
  }
  if(!(action %in% actions)){
    stop("Argument 'action' must be either 'include' or 'exclude'")
  }
  
  #Instantiate column header of interest
  term <- paste(disease,CNV,filt,"perm_p",sep=".")
  if(disease == "ANY_DISEASE" | CNV == "ANY_CNV" | filt == "either_filter"){
    term <- paste("MIN",term,sep=".")
  }
  
  #Get corresponding column index
  cidx <- which(names(df)==term)
  if(length(cidx)==0){
    stop("No column headers match your specified criteria. Please check your input file.")
  }
  
  #Slice df
  if(action=="include"){
    results <- df[which(df[,cidx]<=0.05),]
  }else{
    results <- df[which(df[,cidx]>0.05 | is.na(df[,cidx])),]
  }
  
  #Return results
  names(results)[1] <- "#chr"
  return(results)
}

#list of Rscript options
option_list <- list(
  make_option(c("-C", "--CNV"), type="character", default="ANY_CNV",
              help="CNV class [default '%default']", 
              metavar="character"),
  make_option(c("-d", "--disease"), type="character", default="ANY_DISEASE",
              help="disease cohort [default '%default']", 
              metavar="character"),
  make_option(c("-f", "--filter"), type="character", default="either_filter",
              help="coding+noncoding or noncoding-only [default '%default']", 
              metavar="character"),
  make_option(c("-a", "--action"), type="character", default="include",
              help="include or exclude bins based on criteria [default '%default']", 
              metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="/dev/stdout",
              help="write output to file [default stdout]", 
              metavar="character")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] burden_file",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#Checks for appropriate positional arguments
if(length(args$args) != 1) {
  print_help(parser)
  stop("Must specify an input burden file")
}

#Loads data frame
df <- read.table(args[1],comment.char="",header=T)

#Filters
output <- filter.table(df,
                       disease=opts$disease,
                       CNV=opts$CNV,
                       filt=opts$filt,
                       action=opts$action)
names(output)[1] <- "#chr"

#Write out
write.table(output,opts$outfile,col.names=T,row.names=F,sep="\t",quote=F)
