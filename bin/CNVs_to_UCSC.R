#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Script to send relevant CNVs as a track to a UCSC browser session

##Command line arguments:
# CNVDIR: full path to directory containing CNVs
##Command line options:
# -p/--pheno: phenotype (see external docs for abbreviations)
# -c/--CNV: CNV class (CNV/DEL/DUP)
# -v/--VF: variant frequency (E2/E3/E4/N1)
# -f/--filt: context filter (all/coding/haplosufficient/noncoding/intergenic)
# -g/--gene: query based on gene symbol
# -G/--Gencode: path to gencode gene coordinate BED file (required for gene-based query)
# -l/--locus: query based on locus coordinates
# -b/--buffer: pad query window by ±buffer

#Note: must supply either gene-based query or locus-based query, but not both

#Dev parameters for ERISOne
CNVDIR <- "/data/talkowski/Samples/rCNVmap/data/CNV/CNV_MASTER/"
pheno <- "NDD"
CNV <- "DUP"
VF <- "E4"
filt <- "all"
gene <- "EMB"
geneCoords <- "/data/talkowski/Samples/rCNVmap/data/master_annotations/gencode/gencode.v19.gene_boundaries.protein_coding.bed"
buffer <- 50000

###################
#####Set parameters
###################
options(scipen=1000,stringsAsFactors=F)

######################
#####Load requirements
######################
require("optparse")
require("rtracklayer")

##########################################
#####Helper function to read & filter CNVs
##########################################
importCNVs <- function(pheno,CNV,VF,filt,chr,start,end){
  #Import CNVs
  CNVs <- read.table(paste(CNVDIR,"/",pheno,"/",pheno,".",CNV,".",VF,".GRCh37.",filt,".bed.gz",sep=""))
  names(CNVs)[1:3] <- c("chr","start","end")
  CNVs$chr <- gsub("chr","",CNVs$chr)

  #Filter CNVs
  CNVs <- CNVs[which(as.character(CNVs$chr) == as.character(chr) &
                       as.numeric(CNVs$start) <= as.numeric(end) &
                       as.numeric(CNVs$end) >= as.numeric(start)),]

  #Return CNVs
  return(CNVs)
}

#################################################################
#####Helper function to convert CNV BED file to UCSC track format
#################################################################
formatCNVs <- function(CNVs,color){
  #Format CNV BED
  CNVs.bed <- CNVs[,1:3]
  CNVs.bed[,1] <- paste("chr",CNVs.bed[,1],sep="")
  CNVs.bed[,4] <- paste(round((CNVs[,3]-CNVs[,2])/1000,digits=1),"kb",sep="")
  CNVs.bed[,5] <- 0
  CNVs.bed[,6] <- "+"
  CNVs.bed[,7:8] <- CNVs.bed[2:3]
  CNVs.bed[,9] <- color
  colnames(CNVs.bed) <- c("chrom","chromStart","chromEnd",
                          "name","score","strand",
                          "thickStart","thickEnd","itemRgb")

  #Return formatted CNV bed
  return(CNVs.bed)
}

####################################################################
#####Helper function to generate UCSC view from two sets of CNV BEDs
####################################################################
viewUCSC <- function(CNV.all){
  #Creates CNV track
  CNVRanges <- IRanges(CNV.all$chromStart,CNV.all$chromEnd)
  CNVTrack <- with(CNV.all,
                   GRangesForUCSCGenome("hg19",chrom,CNVRanges,strand,name,
                                        score,thickStart,thickEnd,itemRgb))

  #Starts new browser session
  session <- browserSession("UCSC")

  #Lays track
  track(session,"CNVs",itemRgb=TRUE) <- CNVTrack

  #Sets window
  view.window.df <- data.frame(paste("chr",chr,sep=""),start,end)
  colnames(view.window.df) <- c("chr","start","end")
  view.window <- makeGRangesFromDataFrame(view.window.df,ignore.strand=T)

  #Generates UCSC view
  view <- browserView(session,view.window,pack="CNVs")

  #Returns UCSC view
  return(view)
}

#######################################
#####Rscript command line functionality
#######################################
#List of Rscript options
option_list <- list(
  make_option(c("-p", "--pheno"), type="character", default="GERM",
              help="disease phenotype [default: GERM]",
              metavar="character"),
  make_option(c("-c", "--CNV"), type="character", default="CNV",
              help="CNV class [default: CNV]",
              metavar="character"),
  make_option(c("-v", "--VF"), type="character", default="E2",
              help="variant frequency [default: E2]",
              metavar="character"),
  make_option(c("-f", "--filter"), type="character", default="all",
              help="genomic context filter [default: all]",
              metavar="character"),
  make_option(c("-g", "--gene"), type="character", default=NA,
              help="query based on gene symbol [default: NA]",
              metavar="character"),
  make_option(c("-G", "--Gencode"), type="character", default=NA,
              help="path to Gencode gene coordinate BED [default: NA]",
              metavar="character"),
  make_option(c("-l", "--locus"), type="character", default=NA,
              help="query based on coordinates [default: NA]",
              metavar="character"),
  make_option(c("-b", "--buffer"), type="integer", default=50000,
              help="buffer query coordinates with additional bases [default: 50,000bp]",
              metavar="character")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] CNVDIR",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#Assigns positional arguments to named variables
CNVDIR <- args$args[1]

#Checks for appropriate positional arguments
if(length(args$args) != 1) {
  print_help(parser)
  stop("Must provide path to directory with CNV files")
}

#Assigns options to named variables
pheno <- opts$pheno
CNV <- opts$CNV
VF <- opts$VF
filt <- opts$filter
gene <- opts$gene
geneCoords <- opts$gencode
locus <- opts$locus
buffer <- opts$buffer

#Checks for either gene-based or locus-based query
if(is.na(gene) & is.na(locus)){
  stop("Must supply either query gene [-g] or query locus [-l]")
}
if(!is.na(gene) & !is.na(locus)){
  stop("Must supply either query gene [-g] or query locus [-l], but not both")
}

#Splits coordinates from locus-based query
if(!is.na(locus)){
  chr <- strsplit(locus,split=":")[[1]][1]
  chr <- gsub("chr","",chr)
  start <- as.numeric(strsplit(strsplit(locus,split=":")[[1]][2],"-")[[1]][1])
  end <- as.numeric(strsplit(strsplit(locus,split=":")[[1]][2],"-")[[1]][2])
}

#Gathers coordinates from gene-based query
if(!is.na(gene)){
  geneCoords.file <- read.table(geneCoords,header=F)
  names(geneCoords.file) <- c("chr","start","end","gene")
  chr <- geneCoords.file[which(geneCoords.file$gene==gene),1]
  chr <- gsub("chr","",chr)
  start <- geneCoords.file[which(geneCoords.file$gene==gene),2]
  end <- geneCoords.file[which(geneCoords.file$gene==gene),3]
}

#Pads coordinates with buffer
window.start <- max(start-buffer,0)
window.end <- end+buffer

#Import & filter CNVs
CNV.case <- importCNVs(pheno,CNV,VF,filt,chr,start,end)
CNV.control <- importCNVs("CTRL",CNV,VF,filt,chr,start,end)

#Format CNVs as single UCSC track
CNV.case <- formatCNVs(CNV.case,"red")
CNV.control <- formatCNVs(CNV.control,"gray20")
CNV.all <- rbind(CNV.case,CNV.control)

#Send to UCSC
viewUCSC(CNV.all)

