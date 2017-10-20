#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate locus plot for SMARCA2

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
chr <- "chr9"
start <- 1260000
end <- 2900000
distal.start <- 80050742
distal.end <- 80197059
proximal.start <- 2650652
proximal.end <- 2794964

########################
#####Prepare GViz tracks
########################

#Basic genome coordinates & ideogram
axisTrack <-  GenomeAxisTrack(fontcolor="black",distFromAxis=1.25,labelPos="above",
                              range=IRanges(start=c(distal.start,proximal.start),
                                            end=c(distal.end,proximal.end),
                                            names=c("Distal Region","Proximal Region"),
                              showId=T,col.id="white",col.range=NA,fill.range="red")
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
CNVtracks <- lapply(list("all"),function(filt){
  lapply(list("DEL"),function(CNV){
    dat <- read.table(paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                            CNV,"_",filt,"_density.5kb_bins.bed",sep=""),
                      header=T,comment.char="")
    names(dat)[1] <- "chr"
    gr <- with(dat, GRanges(chr, IRanges(start, end), "CTRL"=CTRL, "NDD"=NDD))
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

#Load fetal brain expression
expression.dat <- read.table(paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                                   "fetal_brain_expression.1kb_bins.bed",sep=""),
                             header=T,comment.char="")
names(expression.dat)[1] <- "chr"
expression.gr <- with(expression.dat, GRanges(chr, IRanges(start, end), "Expression"=Expression))
expression.track <-  DataTrack(range=expression.gr,type="polygon",name="Fetal Brain Expression",
                               col.mountain="black",fill.mountain="black",fill="black")

#Load annotated enhancers
annotated.enhancers <- bed_to_granges(paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                                            "enhancer_elements.bed",sep=""))
annotated.enhancer.track <- AnnotationTrack(annotated.enhancers,name="Fetal Brain Enhancers",
                                            stacking="dense",col.line=cols.NEURO[1],fill=cols.NEURO[1])

#Load fetal brain DHS
DHS.track <- DataTrack(range=paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                                   "fetal_brain_DHS.bedGraph",sep=""),type="polygon",
                       fill=cols.NEURO[1],name="Fetal Brain DHS")

#Load adult brain H3K27ac ChIP-seq track
H3K27ac.track <- DataTrack(range=paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                                       "adult_brain_H3K27ac.bedGraph",sep=""),type="polygon",
                           fill=cols.NEURO[1],name="Adult Brain H3K27ac")

#Load fetal brain H3K4me1 ChIP-seq track
H3K4me1.track <- DataTrack(range=paste(WRKDIR,"plot_data/ExampleLocusPlots/SMARCA2/SMARCA2_locus.",
                                       "fetal_brain_H3K4me1.bedGraph",sep=""),type="polygon",
                           fill=cols.NEURO[1],name="Fetal Brain H3K4me1")



##############
#####Plot view
##############
#Prepare plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/ExampleLoci/SMARCA2/SMARCA2_locus.master.pdf",sep=""),
    height=4,width=8)
#Plot tracks
plotTracks(from=start,to=end,
           c(ideoTrack,axisTrack,unlist(CNVtracks),grTrack,expression.track,
             annotated.enhancer.track,DHS.track,H3K27ac.track,H3K4me1.track),
           background.title="white",col.title="black",col.line="black",col.axis="black",cex.axis=0.5)
#Close device
dev.off()








