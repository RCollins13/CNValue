#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmaps of rCNV signif gene jaccard indexes between pheno groups

###################
#####Set parameters
###################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")
phenos.filenames <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCZ","BEHAV","SEIZ",
                      "HYPO","UNK","SOMA","DRU","GRO","INT","SKEL","HEAD","MSC","CARD",
                      "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSK",
                      "CGEN","CHNK","CBST","CLNG","CEND")
phenos.reorder <- c(1,12,
                    2,3,5,7,8,4,10,11,9,6,
                    13,18,15,20,17,14,19,21,16,22,
                    23,31,28,27,29,25,34,33,35,32,24,30,26)
options(stringsAsFactors=F)

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

###################
#####Load libraries
###################
require(plotrix)

##################################################################
#####Helper function to compute jaccard indexes between two phenos
##################################################################
jaccard <- function(phenoA,phenoB,CNV="BOTH"){
  #Read gene lists for phenoA & phenoB
  if(CNV=="DEL"){
    genesA <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoA,
                                         "_DEL_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                   header=F)[,1])
    genesB <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoB,
                                         "_DEL_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                   header=F)[,1])
  }else if(CNV=="DUP"){
    genesA <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoA,
                                         "_DUP_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                   header=F)[,1])
    genesB <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoB,
                                         "_DUP_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                   header=F)[,1])
  }else{
    genesA.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoA,
                                             "_DEL_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                       header=F)[,1])
    genesB.DEL <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoB,
                                             "_DEL_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                       header=F)[,1])
    genesA.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoA,
                                             "_DUP_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                       header=F)[,1])
    genesB.DUP <- as.vector(read.table(paste(WRKDIR,"plot_data/signif_genes_unique/",phenoB,
                                             "_DUP_E4_exonic.geneScore_Bonferroni_sig.unique.genes.list",sep=""),
                                       header=F)[,1])
    #Merge
    genesA <- unique(union(genesA.DEL,genesA.DUP))
    genesB <- unique(union(genesB.DEL,genesB.DUP))
  }

  #Calculate jaccard index
  jaccard <- length(intersect(genesA,genesB))/length(unique(union(genesA,genesB)))
  return(jaccard)
}

###########################################
#####Compute jaccard matrix for all phenos
###########################################
#Compute jaccard indexes
jacMat <- sapply(phenos.filenames,function(phenoA){
  sapply(phenos.filenames,function(phenoB){
    jaccard(phenoA,phenoB,CNV="BOTH")
  })
})
#Hierarchical clustering of rows/columns
clusters <- hclust(dist(jacMat))
order <- rev(clusters$order)
#Reorder matrix
jacMat <- jacMat[order,order]
#Rename matrix
rownames(jacMat) <- phenos[order]
colnames(jacMat) <- phenos[order]

#########################################
#####Helper function to plot half-heatmap
#########################################
halfHeat <- function(mat,half=T,valMax=0.4,
                     xcolors=NULL,ycolors=NULL,
                     marwd=0.025){
  #Convert df to matrix
  mat <- as.matrix(mat)

  #Round all values down to valMax
  mat[which(mat>valMax)] <- valMax

  #Instatiate margin colors if NULL
  if(is.null(xcolors)){
    xcolors <- rep("white",ncol(mat))
  }
  if(is.null(ycolors)){
    ycolors <- rep("white",nrow(mat))
  }

  #Prepare plotting area
  par(mar=c(3.5,3.5,0.5,0.5),bty="n")
  plot(x=c(-(marwd*ncol(mat)),1.005*ncol(mat)),
       y=c(0.005*ncol(mat),-((marwd+1)*ncol(mat))),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Instantiate color palette
  cols <- rev(rainbow(200)[35:135])
  # cols <- rev(rainbow(150)[1:101])

  #Iterate over cells and plot heatmap
  sapply(1:nrow(mat),function(row){
    sapply(1:ncol(mat),function(col){
      if(half==T & col>=row){
        color <- NA
      }else{
        color <- cols[ceiling(100*(1/valMax)*mat[row,col])+1]
      }
      rect(xleft=col-1,xright=col,
           ybottom=-row,ytop=-row+1,
           border=color,col=color)
    })
  })

  #Add gridlines
  abline(v=1:ncol(mat),h=-1:-nrow(mat),col="white",lwd=0.25)

  #Add outline
  segments(x0=c(0,0),x1=c(nrow(mat)-1,0),
           y0=c(-nrow(mat),-nrow(mat)),
           y1=c(-nrow(mat),-1))
  sapply(1:nrow(mat),function(row){
    sapply(1:ncol(mat),function(col){
      if(half==T & col+1==row){
        segments(x0=c(col-1,col),x1=c(col,col),
                 y0=c(-row+1,-row+1),y1=c(-row+1,-row))
      }
    })
  })

  #Add x-axis margin boxes
  rect(xleft=0:(ncol(mat)-2),
       xright=1:(ncol(mat)-1),
       ytop=-((1+(0.15*marwd))*nrow(mat)),
       ybottom=par("usr")[3],
       col=xcolors,lwd=0.75)

  #Add x-axis margin labels
  sapply(1:(ncol(mat)-1),function(i){
    axis(1,at=i-0.5,tick=F,
         labels=colnames(mat)[i],las=2,line=-0.9)
  })

  #Add y-axis margin boxes
  rect(ybottom=-2:-nrow(mat),
       ytop=-1:-(nrow(mat)-1),
       xright=-0.15*marwd*ncol(mat),
       xleft=par("usr")[1],
       col=xcolors[-1],lwd=0.75)

  #Add y-axis margin labels
  sapply(2:nrow(mat),function(i){
    axis(2,at=-i+0.5,tick=F,
         labels=rownames(mat)[i],las=2,line=-0.9)
  })
}

##########################
#####Generate half heatmap
##########################
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/signifGenes_phenoVsPheno_jaccard.heatmap.pdf",sep=""),
    width=7,height=7)
halfHeat(jacMat,
         xcolors=cols.allPhenos[order],
         ycolors=cols.allPhenos[order])
dev.off()

#######################
#####Plot Jaccard scale
#######################
#Scale for p-value heatmap
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/jaccard_index_blue_green_yellow.scale.pdf",sep=""),
    width=4,height=0.4)
par(mar=c(0.5,0.5,0.5,0.5),bty="o")
plot(x=c(1,101),y=c(0,1),type="n",
     xaxt="n",xlab="",xaxs="i",
     yaxt="n",ylab="",yaxs="i")
cols <- rev(rainbow(200)[35:135])
rect(xleft=0:100,xright=1:101,
     ybottom=0,ytop=1,
     border=cols,col=cols)
axis(3,at=seq(1,101,25),labels=NA)
dev.off()




