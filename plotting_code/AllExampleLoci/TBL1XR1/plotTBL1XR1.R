#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate locus plot for TBL1XR1/NLGN1

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(Gviz)

#####Helper function to read BED data into a GRanges object####
#Note: credit (https://davetang.org/muse/2015/02/04/bed-granges/)
bed_to_granges <- function(file){
  require(GenomicRanges)
  df <- read.table(file,header=F,stringsAsFactors=F)

  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }

  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

#####Set variables for plotting
ref <- "hg19"
chr <- "chr3"
start <- 172579763
end <- 178680636
hot.left.start <- 174145797
hot.left.end <- 175072055
hot.center.start <- 175759747
hot.center.end <- 176688171
hot.right.start <- 176928607
hot.right.end <- 177967562
all.interval.genes <- c("SPATA16","NLGN1","NAALADL2","TBL1XR1","KCNMB2")
relevant.genes <- c("NLGN1","TBL1XR1")

########################
#####Prepare GViz tracks
########################

#Basic genome coordinates & ideogram
axisTrack <-  GenomeAxisTrack(fontcolor="black",distFromAxis=1.25,labelPos="above",
                              range=IRanges(start=c(hot.left.start,hot.center.start,hot.right.start),
                                            end=c(hot.left.end,hot.center.end,hot.right.end),
                                            names=c("Proximal","Medial","Distal")),
                              showId=T,col.id="white",col.range=NA,
                              fill.range=c("red","darkorchid4","blue"))
ideoTrack <- IdeogramTrack(genome=ref,chromosome=chr,
                           col=NA,fill=adjustcolor(cols.NEURO[1],alpha=0.5))

#Prepare genes track
grTrack <- BiomartGeneRegionTrack(genome=ref,chr=chr,start=start,end=end,
                                  collapseTranscripts="meta",transcriptAnnotation="symbol",
                                  col.line="black",col="black",fill="black",
                                  filter=list(with_hgnc=T,with_protein_id=T),
                                  name="Genes")

#Get gene symbols in region
genes <- unique(symbol(grTrack))

#Load CNVs & convert to heatmap data tracks
#Round all CNV frequencies to <0.1%
CNVtracks <- lapply(list("coding","noncoding"),function(filt){
  lapply(list("DEL","DUP"),function(CNV){
    dat <- read.table(paste(WRKDIR,"plot_data/ExampleLocusPlots/TBL1XR1/TBL1XR1_locus.",
                            CNV,"_",filt,"_density.5kb_bins.bed",sep=""),
                      header=T,comment.char="")
    names(dat)[1] <- "chr"
    dat[,4:6] <- apply(dat[,4:6],2,function(vals){
      vals[which(vals>10)] <- 10
      return(vals)
    })
    gr <- with(dat, GRanges(chr, IRanges(start, end), "CTRL"=CTRL, "GERM"=GERM, "CNCR"=CNCR))
    if(CNV=="DEL"){
      gradCol <- c("white","red")
    }else{
      gradCol <- c("white","blue")
    }
      return(DataTrack(gr,ncolor=100,type="heatmap",gradient=gradCol,
             name=paste(CNV," (",filt,")",sep=""),
             showSampleNames=T,col.sampleNames="black"))
  })
})


##############
#####Plot view
##############
#Prepare plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/ExampleLoci/TBL1XR1/TBL1XR1_NLGN1_locus.master.pdf",sep=""),
    height=4,width=8)
#Plot tracks
plotTracks(from=start,to=end,
           c(ideoTrack,axisTrack,grTrack,unlist(CNVtracks)),
           background.title="white",col.title="black",col.line="black",col.axis="black",cex.axis=0.5)
#Close device
dev.off()








