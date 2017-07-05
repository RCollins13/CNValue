#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmap and dendograms of gene sets

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(flashClust)

#####Read correlation matrix
mat <- read.table(paste(WRKDIR,"plot_data/figure3/gene_set_overlaps.matrix.txt",sep=""),header=T)
gsets <- mat[,1]
mat <- mat[,-1]

#####Read geneset groupings

####################################
#####Master function to plot heatmap
####################################
plotHeat <- function(mat,xcols,ycols,buffer=10,squareSize=15){
  #Prepare plotting area
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(x=c(-(buffer+squareSize),ncol(mat)),
       y=c(buffer+squareSize,-(nrow(mat)-1)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Instantiate color palette
  cols <- rev(rainbow(200)[35:135])
  # cols <- rev(rainbow(150)[1:101])

  #Iterate over cells and plot
  sapply(1:nrow(mat),function(row){
    sapply(1:ncol(mat),function(col){
      rect(xleft=row-1,xright=row,
           ybottom=-col,ytop=-col+1,
           border=NA,col=cols[ceiling(100*mat[row,col])+1])
    })
  })

  #Add squares for X legend
  rect(xleft=0:(ncol(mat)-1),xright=1:ncol(mat),
       ybottom=buffer,ytop=buffer+squareSize,
       border=NA,col=xcols)

  #Add squares for Y legend
  rect(xleft=-(buffer+squareSize),xright=-buffer,
       ybottom=-ncol(mat):-1,ytop=-(ncol(mat)-1):0,
       border=NA,col=ycols)
}

########################
#####Unclustered heatmap
########################
# plotHeat(mat)

########################
#####Clustered heatmap
########################
xclusters <- hclust(dist(t(mat)))
plot(xclusters)
xgroups <- cutree(xclusters,6)
yclusters <- hclust(dist(mat))
plot(yclusters)
ygroups <- cutree(yclusters,6)
mat.ordered <- mat[xclusters$order,yclusters$order]
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/gene_set_overlap.png",sep=""),
    width=4,height=4,units="in",res=1000)
plotHeat(mat.ordered)
dev.off()



