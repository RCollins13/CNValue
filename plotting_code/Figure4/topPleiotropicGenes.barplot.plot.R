#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate dotplot grid of top 25 most pleiotropic genes

###################
#####Set parameters
###################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")
phenos.reorder <- c(1,12,2,3,5,7,8,4,10,11,9,6,
                    13,18,15,20,17,14,19,21,16,22,
                    23,31,28,27,29,25,34,33,35,32,
                    24,30,26)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")
cols.allPhenos <- c(cols.GERM[1],
                    rep(cols.NEURO[1],10),
                    cols.GERM[1],
                    rep(cols.SOMA[1],10),
                    rep(cols.CNCR[1],13))

####################################
#####Helper function to plot heatmap
####################################
plotHeat <- function(mat,
                     xcolors=NULL,
                     ylabs=NULL,
                     marwd=0.001){
  #Convert df to matrix
  mat <- as.matrix(mat)

  #Instatiate margin colors if NULL
  if(is.null(xcolors)){
    xcolors <- rep("white",ncol(mat))
  }

  #Prepare plotting area
  par(mar=c(0.5,3.5,3.5,0.5),bty="n")
  plot(x=c(0,ncol(mat)),y=c(-nrow(mat),marwd*ncol(mat)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Iterate over cells and plot heatmap
  sapply(1:nrow(mat),function(row){
    sapply(1:ncol(mat),function(col){
      if(mat[row,col]=="DEL"){
        color <- "red"
      }else if(mat[row,col]=="DUP"){
        color <- "blue"
      }else if(mat[row,col]=="BOTH"){
        color <- cols.CTRL[1]
      }else{
        color <- NA
      }
      rect(xleft=col-1,xright=col,
           ybottom=-row,ytop=-row+1,
           border=NA,col=color)
    })
  })

  #Add gridlines
  abline(v=1:ncol(mat),h=-1:-nrow(mat),col="white",lwd=0.25)

  #Add x-axis margin boxes
  rect(xleft=0:(ncol(mat)-2),
       xright=1:(ncol(mat)-1),
       ybottom=0,ytop=par("usr")[4],
       col=xcolors,lwd=0.75)

  #Add x-axis margin labels
  sapply(1:(ncol(mat)-1),function(i){
    axis(3,at=i-0.5,tick=F,cex.axis=1.1,
         labels=colnames(mat)[i],las=2,line=-0.9)
  })

  #Add y-axis margin labels
  if(!is.null(ylabs)){
    sapply(2:nrow(mat),function(i){
      axis(2,at=-i+0.5,tick=F,cex.axis=1.1,
           labels=ylabs[i],las=2,line=-0.9)
    })
  }
}

##############
#####Read data
##############
df <- read.table(paste(WRKDIR,"plot_data/figure4/topGenes_assocByPheno.txt",sep=""),header=T)
genes <- df[,1]
#Count DEL-DUP in GERM and DEL-DUP in CNCR
dGERM <- apply(df[,2:23],1,function(vals){
  DEL <- length(which(vals=="BOTH" | vals=="DEL"))
  DUP <- length(which(vals=="BOTH" | vals=="DUP"))
  d <- (DEL/22)-(DUP/22)
  return(d)
})
dCNCR <- apply(df[,24:ncol(df)],1,function(vals){
  DEL <- length(which(vals=="BOTH" | vals=="DEL"))
  DUP <- length(which(vals=="BOTH" | vals=="DUP"))
  d <- (DEL/13)-(DUP/13)
  return(d)
})
dCross <- dGERM-dCNCR
gene.order <- c(order(dCross,decreasing=T)[1:13],14,
                14+order(dCross[15:18],decreasing=T),19,
                19+order(dCross[20:length(dCross)],decreasing=T))
df <- df[gene.order,
         phenos.reorder+1]

#################
#####Plot heatmap
#################
plotHeat(df,xcolors=cols.allPhenos[phenos.reorder],
         ylabs=genes[gene.order])
