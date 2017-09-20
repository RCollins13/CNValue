#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmaps of rCNV signif gene enrichments in 263 gene sets

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")
phenos.reorder <- c(1,12,
                    2,3,5,7,8,4,10,11,9,6,
                    13,18,15,20,17,14,19,21,16,22,
                    23,31,28,27,29,25,34,33,35,32,24,30,26)

#####Set color vectors
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

#####Load libraries
require(plotrix)

############################################
#####Helper function to compute Fisher stats
############################################
fisherStats <- function(df){
  #Run Fisher tests
  fisher.res <- as.data.frame(t(sapply(1:nrow(df),function(r){
    #Subset row
    vals <- as.numeric(df[r,-1])

    #Iterate over all phenotypes
    fisher.res <- as.vector(sapply(seq(3,length(vals),2),function(i){
      fisher.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                       vals[i]-vals[i+1],vals[i+1])),nrow=2)
      fisher.res <- fisher.test(x=fisher.df,alternative="greater")
      return(as.numeric(c(fisher.res$estimate,fisher.res$conf.int,fisher.res$p.value)))
    }))

    #Return values
    return(fisher.res)
  })))

  #Assign names to results
  # rownames(fisher.res) <- df[,1]
  names(fisher.res) <- sapply(phenos,function(pheno){
    return(paste(pheno,c(".OR",".lowerCI",".upperCI",".p"),sep=""))
  })

  #Return results
  return(fisher.res)
}

##############################################
#####Helper function to compute binomial stats
##############################################
binomStats <- function(df){
  #Run binomial tests
  binom.res <- as.data.frame(t(sapply(1:nrow(df),function(r){
    #Subset row
    vals <- as.numeric(df[r,-1])

    #Iterate over all phenotypes
    binom.res <- as.vector(sapply(seq(3,length(vals),2),function(i){
      res <- binom.test(x=vals[i+1],n=vals[i],p=vals[2]/vals[1])
      return(as.numeric(c(res$estimate/as.numeric(res$null.value),
                          res$conf.int[1]/as.numeric(res$null.value),
                          res$conf.int[2]/as.numeric(res$null.value),
                          res$p.value)))
    }))

    #Return values
    return(binom.res)
  })))

  #Assign names to results
  # rownames(binom.res) <- df[,1]
  names(binom.res) <- sapply(phenos,function(pheno){
    return(paste(pheno,c(".fold",".lowerCI",".upperCI",".p"),sep=""))
  })

  #Return results
  return(binom.res)
}

##############################################################
#####Helper function to read data & compute all relevant stats
##############################################################
readData <- function(CNV,VF,context,sig,constrained=T){
  #Read file
  if(constrained==T){
    df <- read.table(paste(WRKDIR,"plot_data/signif_genes_geneset_comparisons/allPhenos_allSets/",
                           CNV,"_",VF,"_",context,"_",sig,".comparisons.txt",sep=""),
                     header=T,comment.char="")
  }else{
    df <- read.table(paste(WRKDIR,"plot_data/signif_genes_geneset_comparisons/allPhenos_allSets_noConstrained/",
                           CNV,"_",VF,"_",context,"_",sig,".comparisons.txt",sep=""),
                     header=T,comment.char="")
  }
  colnames(df)[1] <- "geneset"

  #Compute Fisher stats
  fisher.res <- fisherStats(df)

  #Compute binomial p-values
  binom.res <- binomStats(df)

  #Prepare output list
  results <- list("counts"=df,
                  "fisher"=fisher.res,
                  "binom"=binom.res)
  return(results)
}

###############################################################
#####Helper function to extract p-values from a fisher/binom df
###############################################################
extractP <- function(df,bonf=F,log=T){
  #Assume p-values are stored in every fourth column
  df.p <- df[,seq(4,ncol(df),4)]

  #-log10-transform p-values
  df.p <- -log10(df.p)

  #Add headers
  names(df.p) <- phenos

  #Return
  return(df.p)
}

#########################################
#####Helper function to compute OR deltas
#########################################
deltaOR <- function(dfA,dfB){
  res <- sapply(1:nrow(dfA),function(row){
    sapply(1:ncol(dfA),function(col){
      if(!is.infinite(dfA[row,col]) & !is.infinite(dfB[row,col])){
        return(dfA[row,col]-dfB[row,col])
      }else{
        return(0)
      }
    })
  })

  return(t(res))
}

##########################################################
#####Helper function to extract ORs from a fisher/binom df
##########################################################
extractOR <- function(df,bonf=F,log=T){
  #Assume ORs are stored in every first column
  df.OR <- df[,seq(1,ncol(df),4)]

  #log2-transform p-values
  df.OR <- log2(df.OR)

  #Add headers
  names(df.OR) <- phenos

  #Return
  return(df.OR)
}

################################################
#####Helper function to generate p-value heatmap
################################################
p.heat <- function(df,pmax=10,bonf=F,
                   xcluster=F,ycluster=T,
                   xorder=NULL,yorder=NULL,
                   xcolors=NULL,ycolors=NULL,
                   marwd=0.025,xcats=NULL,
                   lv=NULL,outline=F){
  #Convert df to matrix
  mat <- as.matrix(df)

  #Round down all values to pmax
  mat[which(mat>pmax)] <- pmax

  #Instatiate margin colors if NULL
  if(is.null(xcolors)){
    xcolors <- rep("white",ncol(mat))
  }
  if(is.null(ycolors)){
    ycolors <- rep("white",nrow(mat))
  }

  #Only retain rows with at least one Bonferroni-corrected significant result (if optioned)
  p.cutoff <- -log10(0.05/length(mat))
  if(bonf==T){
    keepers <- apply(mat,1,function(vals){
      return(any(vals>p.cutoff))
    })
    keepers <- which(keepers)
  }else{
    keepers <- 1:nrow(mat)
  }
  mat <- mat[keepers,]
  ycolors <- ycolors[keepers]

  #Cluster x-axis (if optioned)
  if(is.null(xorder)){
    if(xcluster==T){
      xclusters <- hclust(dist(t(mat)))
      xorder <- rev(xclusters$order)
    }else{
      xorder <- 1:ncol(mat)
    }
  }

  #Cluster y-axis (if optioned)
  if(is.null(yorder)){
    if(ycluster==T){
      yclusters <- hclust(dist(mat))
      yorder <- rev(yclusters$order)
    }else{
      yorder <- 1:nrow(mat)
    }
  }

  #Reorder mat & colors by clustering
  mat.reordered <- mat[yorder,xorder]
  xcolors <- xcolors[xorder]
  ycolors <- ycolors[yorder]
  xcats <- xcats[xorder]
  keepers <- keepers[yorder]

  #Prepare plotting area
  par(mar=c(0.5,1.5,0.5,1),bty="n")
  plot(x=c(-(marwd*ncol(mat.reordered)),1.005*ncol(mat.reordered)),
       y=c(((marwd+0.15)*nrow(mat.reordered)),-nrow(mat.reordered)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Instantiate color palette
  # cols <- rev(rainbow(200)[35:135])
  cols <- colorRampPalette(c("#08176b","#138934","#ffff59"))(1001)
  # cols <- rev(rainbow(150)[1:101])

  #Iterate over cells and plot heatmap
  sapply(1:nrow(mat.reordered),function(row){
    sapply(1:ncol(mat.reordered),function(col){
      color <- cols[ceiling(100*mat.reordered[row,col])+1]
      rect(xleft=col-1,xright=col,
           ybottom=-row,ytop=-row+1,
           border=NA,col=color)
    })
  })

  #Add vertical lines
  if(!is.null(lv)){
    abline(v=lv,col="white")
  }

  #Add thin white gridlines
  abline(h=-nrow(mat):nrow(mat),
         v=-nrow(mat):nrow(mat),
         col="white",lwd=0.25)

  #Add x-axis margin boxes
  rect(xleft=0:(ncol(mat.reordered)-1),
       xright=1:ncol(mat.reordered),
       ybottom=0.15*marwd*nrow(mat.reordered),
       ytop=marwd*nrow(mat.reordered),
       col=xcolors,lwd=0.75)

  #Add x-axis margin labels
  sapply(1:ncol(mat.reordered),function(i){
    text(x=i-1,y=1.1*marwd*nrow(mat.reordered),
         labels=xcats[i],srt=90,pos=4,cex=0.75)
  })

  #Add y-axis margin colors
  rect(xleft=-marwd*ncol(mat.reordered),
       xright=-0.15*marwd*ncol(mat.reordered),
       ybottom=-1:-nrow(mat.reordered),
       ytop=0:(-nrow(mat.reordered)+1),
       col=ycolors,lwd=0.75)

  #Add y-axis numbering
  axis(2,at=seq(-5,-nrow(mat.reordered),-5)+0.5,
       tck=-0.01,lwd=0.75,labels=NA)
  axis(2,at=seq(-5,-nrow(mat.reordered),-5)+0.5,
       tick=F,line=-0.6,las=2,cex.axis=0.9,
       labels=seq(5,nrow(mat.reordered),5))

  #Add right y-axis asterisks for categories with >0 Bonf-sig groups
  sig <- apply(mat.reordered,1,function(vals){
    return(any(vals>p.cutoff))
  })
  sig <- which(sig)
  sig.chunks <- 6
  sapply(1:sig.chunks,function(i){
    axis(4,at=-sig[seq(i,length(sig),sig.chunks)]+0.5,
         tick=F,line=-0.9,cex.axis=1.6,font=2,col.axis="red",
         labels=rep("*",length(seq(i,length(sig),sig.chunks))))
  })

  #Iterate over cells and plot red borders (if optioned)
  if(outline==T){
    sapply(1:nrow(mat.reordered),function(row){
      sapply(1:ncol(mat.reordered),function(col){
        if(mat.reordered[row,col]>p.cutoff){
          rect(xleft=col-1,xright=col,
               ybottom=-row,ytop=-row+1,
               border="red",col=NA,lwd=2)
        }
      })
    })
  }

  #Add cleanup rectangle around border
  rect(xleft=0,xright=ncol(mat.reordered),
       ytop=0,ybottom=-nrow(mat.reordered),
       col=NA,border="black")

  #Return column numbers of keepers
  return(keepers)
}

##############################################
#####Helper function to generate delta heatmap
##############################################
delta.heat <- function(df,pmax=3,thresh=NA,
                   xcluster=F,ycluster=T,
                   xorder=NULL,yorder=NULL,
                   xcolors=NULL,ycolors=NULL,
                   marwd=0.025,xcats=NULL,
                   lv=NULL){
  #Convert df to matrix
  mat <- as.matrix(df)

  #Fit all values on ~ [-pmax,pmax]
  mat[which(mat>pmax)] <- pmax
  mat[which(mat< -pmax)] <- -pmax

  #Instatiate margin colors if NULL
  if(is.null(xcolors)){
    xcolors <- rep("white",ncol(mat))
  }
  if(is.null(ycolors)){
    ycolors <- rep("white",nrow(mat))
  }

  #Only retain rows with values above thresh (if optioned)
  if(!is.na(thresh)){
    keepers <- apply(mat,1,function(vals){
      return(any(vals>thresh))
    })
    keepers <- which(keepers)
  }else{
    keepers <- 1:nrow(mat)
  }
  mat <- mat[keepers,]
  ycolors <- ycolors[keepers]

  #Cluster x-axis (if optioned)
  if(is.null(xorder)){
    if(xcluster==T){
      xclusters <- hclust(dist(t(mat)))
      xorder <- rev(xclusters$order)
    }else{
      xorder <- 1:ncol(mat)
    }
  }

  #Cluster y-axis (if optioned)
  if(is.null(yorder)){
    if(ycluster==T){
      yclusters <- hclust(dist(mat))
      yorder <- rev(yclusters$order)
    }else{
      yorder <- 1:nrow(mat)
    }
  }

  #Reorder mat & colors by clustering
  mat.reordered <- mat[yorder,xorder]
  xcolors <- xcolors[xorder]
  ycolors <- ycolors[yorder]
  xcats <- xcats[xorder]
  keepers <- keepers[yorder]

  #Prepare plotting area
  par(mar=c(0.5,1.5,0.5,1),bty="n")
  plot(x=c(-(marwd*ncol(mat.reordered)),1.005*ncol(mat.reordered)),
       y=c(((marwd+0.15)*nrow(mat.reordered)),-nrow(mat.reordered)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Instantiate color palette
  cols <- colorRampPalette(c("blue","white","red"))((20*pmax)+1)
  # cols <- rev(rainbow(150)[1:101])

  #Iterate over cells and plot heatmap
  sapply(1:nrow(mat.reordered),function(row){
    sapply(1:ncol(mat.reordered),function(col){
      color <- cols[ceiling(10*mat.reordered[row,col])+(10*pmax)+1]
      rect(xleft=col-1,xright=col,
           ybottom=-row,ytop=-row+1,
           border=NA,col=color)
    })
  })

  #Add vertical lines
  if(!is.null(lv)){
    abline(v=lv,col="white")
  }

  #Add thin white gridlines
  abline(h=-nrow(mat):nrow(mat),
         v=-nrow(mat):nrow(mat),
         col="white",lwd=0.25)

  #Add x-axis margin boxes
  rect(xleft=0:(ncol(mat.reordered)-1),
       xright=1:ncol(mat.reordered),
       ybottom=0.15*marwd*nrow(mat.reordered),
       ytop=marwd*nrow(mat.reordered),
       col=xcolors,lwd=0.75)

  #Add x-axis margin labels
  sapply(1:ncol(mat.reordered),function(i){
    text(x=i-1,y=1.1*marwd*nrow(mat.reordered),
         labels=xcats[i],srt=90,pos=4,cex=0.75)
  })

  #Add y-axis margin colors
  rect(xleft=-marwd*ncol(mat.reordered),
       xright=-0.15*marwd*ncol(mat.reordered),
       ybottom=-1:-nrow(mat.reordered),
       ytop=0:(-nrow(mat.reordered)+1),
       col=ycolors,lwd=0.75)

  #Add y-axis numbering
  axis(2,at=seq(-5,-nrow(mat.reordered),-5)+0.5,
       tck=-0.01,lwd=0.75,labels=NA)
  axis(2,at=seq(-5,-nrow(mat.reordered),-5)+0.5,
       tick=F,line=-0.6,las=2,cex.axis=0.9,
       labels=seq(5,nrow(mat.reordered),5))

  #Add cleanup rectangle around border
  rect(xleft=0,xright=ncol(mat.reordered),
       ytop=0,ybottom=-nrow(mat.reordered),
       col=NA,border="black")

  #Return column numbers of keepers
  return(keepers)
}

##################################################################
####Generate heatmap of all gene sets vs all phenotypes for figure
##################################################################
#Read data
DEL <- readData("DEL","E4","exonic","Bonferroni",constrained=T)
DEL.p <- extractP(DEL$binom)
DEL.OR <- extractOR(DEL$binom)
DEL.noCon <- readData("DEL","E4","exonic","Bonferroni",constrained=F)
DEL.noCon.p <- extractP(DEL.noCon$binom)
DEL.noCon.OR <- extractOR(DEL.noCon$binom)
DUP <- readData("DUP","E4","exonic","Bonferroni",constrained=T)
DUP.p <- extractP(DUP$binom)
DUP.OR <- extractOR(DUP$binom)
DUP.noCon <- readData("DUP","E4","exonic","Bonferroni",constrained=F)
DUP.noCon.p <- extractP(DUP.noCon$binom)
DUP.noCon.OR <- extractOR(DUP.noCon$binom)
#Select which columns are significant in DEL and/or DUP
DEL.sig <- apply(DEL.p,1,function(vals){
  return(any(vals>-log10(0.05/(nrow(DEL.p)*ncol(DEL.p)))))
})
DUP.sig <- apply(DUP.p,1,function(vals){
  return(any(vals>-log10(0.05/(nrow(DUP.p)*ncol(DUP.p)))))
})
all.sig <- sort(unique(c(which(DEL.sig),which(DUP.sig))))
# Subset dfs to keep columns that are significant in either
DEL.p.eitherSig <- DEL.p[all.sig,]
DEL.OR.eitherSig <- DEL.OR[all.sig,]
DEL.noCon.p.eitherSig <- DEL.noCon.p[all.sig,]
DEL.noCon.OR.eitherSig <- DEL.noCon.OR[all.sig,]
DUP.p.eitherSig <- DUP.p[all.sig,]
DUP.OR.eitherSig <- DUP.OR[all.sig,]
DUP.noCon.p.eitherSig <- DUP.noCon.p[all.sig,]
DUP.noCon.OR.eitherSig <- DUP.noCon.OR[all.sig,]
#Cluster based on joined DEL & DUP to determine y-ordering
joined <- cbind(DEL.p.eitherSig,DUP.p.eitherSig)
rowOrder <- hclust(dist(as.matrix(joined)))$order
#Compute delta p and delta OR for gene sets of interest
delta.p.eitherSig <- DEL.p.eitherSig - DUP.p.eitherSig
delta.noCon.p.eitherSig <- DEL.noCon.p.eitherSig - DUP.noCon.p.eitherSig
delta.OR.eitherSig <- deltaOR(DEL.OR.eitherSig,DUP.OR.eitherSig)
delta.noCon.OR.eitherSig <- deltaOR(DEL.noCon.OR.eitherSig,DUP.noCon.OR.eitherSig)
#Set colors for y-axis categories of interest
catCols <- c("#EAEBEC","#EAEBEC","red","red","red","red","red",
             "red","red","#00BFF4","#7B2AB3","#00BFF4","#00BFF4",
             "#00BFF4","#00BFF4","#EC008D","#EC008D","#EC008D",
             "#FFCB00","#FFCB00","#FFCB00","#FFCB00","#FFCB00",
             "#FFCB00","#FFCB00","#7B2AB3","#FFCB00","#FFCB00",
             "#FFCB00","#FFCB00","#FF6A09","black","#FF6A09",
             "black","#FF6A09","black","#FF6A09","black","black",
             "black","#FF6A09","#FF6A09","black","black","#FF6A09",
             "#FF6A09","#FF6A09","#FF6A09","#FF6A09","#EAEBEC",
             "#EAEBEC","#EAEBEC","#EAEBEC")
#Plot heatmaps
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DEL_Bonf_wConst.heatmap.pdf",sep=""),
    width=5,height=5)
p.heat(DEL.p.eitherSig,xcluster=F,ycluster=F,bonf=F,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xcats=phenos,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DEL_Bonf_noConst.heatmap.pdf",sep=""),
    width=5,height=5)
p.heat(DEL.noCon.p.eitherSig,xcluster=F,ycluster=F,bonf=F,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xcats=NULL,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DUP_Bonf_wConst.heatmap.pdf",sep=""),
    width=5,height=5)
p.heat(DUP.p.eitherSig,xcluster=F,ycluster=F,bonf=F,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xcats=phenos,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DUP_Bonf_noConst.heatmap.pdf",sep=""),
    width=5,height=5)
p.heat(DUP.noCon.p.eitherSig,xcluster=F,ycluster=F,bonf=F,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DELTA_Bonf_wConst.heatmap.pdf",sep=""),
    width=5,height=5)
delta.heat(delta.p.eitherSig,xcluster=F,ycluster=F,pmax=3,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xcats=phenos,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/DELTA_Bonf_noConst.heatmap.pdf",sep=""),
    width=5,height=5)
delta.heat(delta.noCon.p.eitherSig,xcluster=F,ycluster=F,pmax=3,
       xcolors=cols.allPhenos,
       ycolors=catCols,
       xorder=phenos.reorder,
       yorder=rowOrder)
dev.off()

########################
#####Plot keys & legends
########################
#Dots for gene set key
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/heatmap_categories_key.pdf",sep=""),
    width=1/8,height=length(catCols)/8)
par(mar=rep(0.1,4),bty="n")
plot(x=c(0.5,1.5),y=c(0,-length(catCols)),type="n",
     xaxt="n",xlab="",xaxs="i",
     yaxt="n",ylab="",yaxs="i")
points(x=rep(1,length(catCols)),
       y=(-1:-length(catCols))+0.5,
       pch=22,bg=catCols[rowOrder],lwd=0.8)
dev.off()
#Scale for p-value heatmap
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/blue_green_yellow_p_scale.pdf",sep=""),
    width=4,height=0.4)
par(mar=c(0.5,0.5,0.5,0.5),bty="o")
plot(x=c(0,101),y=c(0,1),type="n",
     xaxt="n",xlab="",xaxs="i",
     yaxt="n",ylab="",yaxs="i")
cols <- colorRampPalette(c("#08176b","#138934","#ffff59"))(101)
rect(xleft=0:100,xright=1:101,
     ybottom=0,ytop=1,
     border=cols,col=cols)
abline(v=10*-log10(c(0.05,0.05/(35*nrow(DEL$counts)))),
       lty=c(2,1))
axis(3,at=10*-log10(0.05),labels=NA)
axis(3,at=10*-log10(0.05/(35*nrow(DEL$counts))),labels=NA)
axis(1,at=seq(0.5,100.5,20),labels=NA)
dev.off()
#Scale for delta heatmap
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/red_white_blue_p_scale.pdf",sep=""),
    width=4,height=0.4)
par(mar=c(0.5,0.5,0.5,0.5),bty="o")
plot(x=c(0,601),y=c(0,1),type="n",
     xaxt="n",xlab="",xaxs="i",
     yaxt="n",ylab="",yaxs="i")
cols <- colorRampPalette(c("blue","white","red"))(601)
rect(xleft=0:600,xright=1:601,
     ybottom=0,ytop=1,
     border=cols,col=cols)
axis(1,seq(0.5,600.6,100),labels=NA)
dev.off()


