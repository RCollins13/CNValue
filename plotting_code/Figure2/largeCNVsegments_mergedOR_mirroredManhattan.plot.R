#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate mirrored manhattan plot of ORs per phenotype for del/dup at all significant loci

####################################
#####Set parameters & load libraries
####################################
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
require(beeswarm)

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

###################################################
#####Helper function to read regular manhattan data
###################################################
readData <- function(pheno,VF,filt){
  df <- read.table(paste(WRKDIR,"plot_data/figure2/",pheno,".",
                         VF,".",filt,".manhattan_pvals.bed",sep=""),header=F)
  names(df) <- c("chr","start","end","pDEL","pDUP")
  df$pDEL <- -log10(df$pDEL)
  df$pDUP <- -log10(df$pDUP)
  df$pos <- apply(df[,2:3],1,mean)
  return(df[,c(1,6,4:5)])
}

##############
#####Read data
##############
#P-values
DEL.p <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DEL_pVal.bed",sep=""),header=F)
colnames(DEL.p) <- c("chr","start","end",phenos)
DEL.p <- DEL.p[,c(1:3,phenos.reorder+3)]
DUP.p <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DUP_pVal.bed",sep=""),header=F)
colnames(DUP.p) <- c("chr","start","end",phenos)
DUP.p <- DUP.p[,c(1:3,phenos.reorder+3)]
#Odds ratios
DEL.OR <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DEL_OR.bed",sep=""),header=F)
colnames(DEL.OR) <- c("chr","start","end",phenos)
DEL.OR <- DEL.OR[,c(1:3,phenos.reorder+3)]
DUP.OR <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DUP_OR.bed",sep=""),header=F)
colnames(DUP.OR) <- c("chr","start","end",phenos)
DUP.OR <- DUP.OR[,c(1:3,phenos.reorder+3)]
#Manhattan base
manhat.base <- readData("NEURO","E2","all")

################################################
######Helper function to plot mirrored manhattan
################################################
mirrorManhattan <- function (DEL.OR,DUP.OR,
                             manhat.base,
                             maxOR=256,infOR=512){
  #Gather list of unique contigs
  contigs <- unique(manhat.base[,1])
  contigs <- contigs[which(!(is.na(contigs)))]

  #Create index data frame
  indexes <- as.data.frame(t(sapply(contigs,function(chr){
    return(c(chr,0,max(manhat.base[which(manhat.base[,1]==chr),2])))
  })))
  indexes$sum <- cumsum(indexes[,3])

  #Create new plotting dfs with modified coordinates & formatted values
  DEL.plot <- t(apply(DEL.OR,1,function(vals){
    #Quick fix to values
    vals <- as.vector(vals)
    vals[1] <- gsub("chr","",vals[1])
    vals[-1] <- as.numeric(vals[-1])

    #Get nearest plotting coordinate
    avgCoord <- mean(as.numeric(vals[2:3]))
    deltaCoords <- avgCoord-manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos
    nearestCoordIdx <- head(which(abs(deltaCoords)==min(abs(deltaCoords))),1)
    nearestCoord <- manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos[nearestCoordIdx]
    plottingCoord <- indexes[which(indexes[,1]==as.character(vals[1])),]$sum-indexes[which(indexes[,1]==as.character(vals[1])),3]+nearestCoord

    #Format ORs
    ORs <- as.numeric(as.vector(vals[-c(1:3)]))
    ORs[which(!is.infinite(ORs) & ORs>maxOR)] <- maxOR
    ORs[which(is.infinite(ORs) & ORs>0)] <- infOR
    ORs[which(is.na(ORs))] <- 0
    ORs <- log2(ORs)
    ORs[which(is.infinite(ORs) & ORs<0)] <- 0
    ORs[which(ORs<0)] <- 0

    #Return formatted values for plotting
    return(as.numeric(c(plottingCoord,ORs)))
  }))
  DUP.plot <- t(apply(DUP.OR,1,function(vals){
    #Quick fix to values
    vals <- as.vector(vals)
    vals[1] <- gsub("chr","",vals[1])
    vals[-1] <- as.numeric(vals[-1])

    #Get nearest plotting coordinate
    avgCoord <- mean(as.numeric(vals[2:3]))
    deltaCoords <- avgCoord-manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos
    nearestCoordIdx <- head(which(abs(deltaCoords)==min(abs(deltaCoords))),1)
    nearestCoord <- manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos[nearestCoordIdx]
    plottingCoord <- indexes[which(indexes[,1]==as.character(vals[1])),]$sum-indexes[which(indexes[,1]==as.character(vals[1])),3]+nearestCoord

    #Format ORs
    ORs <- as.numeric(as.vector(vals[-c(1:3)]))
    ORs[which(!is.infinite(ORs) & ORs>maxOR)] <- maxOR
    ORs[which(is.infinite(ORs) & ORs>0)] <- infOR
    ORs[which(is.na(ORs))] <- 0
    ORs <- log2(ORs)
    ORs[which(is.infinite(ORs) & ORs<0)] <- 0
    ORs[which(ORs<0)] <- 0

    #Return formatted values for plotting
    return(as.numeric(c(plottingCoord,ORs)))
  }))

  #Set height of contig boxes
  ymax <- log2(maxOR)
  boxht <- 0.06*ymax

  #Prepare plotting window
  par(mar=c(0.2,2.5,0.2,0.2),bty="n")
  plot(x=c(-0.01*max(indexes[,4]),1.01*max(indexes[,4])),y=c(-1.02*(ymax+boxht),1.02*(ymax+boxht)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Add background shading & gridlines
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=cols.CTRL[4])
  abline(v=indexes[,4],col="white",lwd=2)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=-boxht,ytop=boxht,col="white",border=NA)
  abline(h=c((1:log2(maxOR))+boxht,(-1:-log2(maxOR))-boxht),col=cols.CTRL[3])
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=c(par("usr")[3],ymax+boxht),
       ytop=c(-ymax-boxht,par("usr")[4]),
       border=NA,col="white")
  abline(h=c(ymax+boxht,-ymax-boxht))
  rect(xleft=c(par("usr")[1],max(indexes[,4])),xright=c(0,par("usr")[2]),
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col="white",border=NA)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA,border="white",lwd=2)

  #Add y-axis
  y.at <- 1:ymax
  #Top y-axis
  axis(2,at=y.at+boxht,labels=NA)
  axis(2,at=y.at+boxht,tick=F,line=-0.3,labels=2^y.at,las=2)
  # axis(2,at=ymax+boxht,labels=NA)
  # axis(2,at=ymax+boxht,tick=F,line=-0.3,labels="Inf",las=2)
  #Bottom y-axis
  axis(2,at=-(y.at+boxht),labels=NA)
  axis(2,at=-(y.at+boxht),tick=F,line=-0.3,labels=2^y.at,las=2)
  # axis(2,at=-ymax-boxht,labels=NA)
  # axis(2,at=-ymax-boxht,tick=F,line=-0.3,labels="Inf",las=2)

  #Determine significant loci
  del.sig <- apply(DEL.p[,-c(1:3)],1,function(p){
    if(any(p<=10^-8)){
      return(1)
    }else{
      return(0)
    }
  })
  del.sig <- which(del.sig>0)
  dup.sig <- apply(DUP.p[,-c(1:3)],1,function(p){
    if(any(p<=10^-8)){
      return(1)
    }else{
      return(0)
    }
  })
  dup.sig <- which(dup.sig>0)

  #Plots points
  apply(DEL.plot[del.sig,],1,function(vals){
    points(x=rep(vals[1],times=length(vals)-1),y=vals[-1]+boxht,
           pch=21,bg=cols.allPhenos[phenos.reorder],cex=0.5,lwd=0)
  })
  apply(DUP.plot[dup.sig,],1,function(vals){
    points(x=rep(vals[1],times=length(vals)-1),y=-vals[-1]-boxht,
           pch=21,bg=cols.allPhenos[phenos.reorder],cex=0.5,lwd=0)
  })
  apply(DEL.plot[del.sig,],1,function(vals){
    medians <- c(mean(vals[2:23]),mean(vals[4:13]),
                 mean(vals[14:23]),mean(vals[24:36]))
    points(x=rep(vals[1],times=4),y=medians+boxht,
           pch=23,bg=c(cols.GERM[1],cols.NEURO[1],
                       cols.SOMA[1],cols.CNCR[1]),
           cex=0.85,lwd=0.5)
  })
  apply(DUP.plot[dup.sig,],1,function(vals){
    medians <- c(mean(vals[2:23]),mean(vals[4:13]),
                 mean(vals[14:23]),mean(vals[24:36]))
    points(x=rep(vals[1],times=4),y=-medians-boxht,
           pch=23,bg=c(cols.GERM[1],cols.NEURO[1],
                       cols.SOMA[1],cols.CNCR[1]),
           cex=0.85,lwd=0.5)
  })

  #Adds chromosome labels
  sapply(1:length(contigs),function(i){
    #Rectangle
    rect(xleft=indexes[i,4]-indexes[i,3],xright=indexes[i,4],
         ybottom=-boxht,ytop=boxht,col="white")
    if(i<=13){
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.85)
    }else if(i<=18){
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.6,srt=90)
    }else{
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.45,srt=90)
    }
  })

  #Add ticks at significant loci
  segments(x0=DEL.plot[del.sig,1],x1=DEL.plot[del.sig,1],
           y0=0.75*boxht,y1=boxht,col="red")
  segments(x0=DUP.plot[dup.sig,1],x1=DUP.plot[dup.sig,1],
           y0=-0.75*boxht,y1=-boxht,col="blue")
}

###########################################
#####Plot median ORs for all sites by pheno
###########################################
#Note: mostly copied from mirrored manhattan code for convenience
plotMedianORs <- function(OR,p,manhat.base,
                          maxOR=256,infOR=512){
  #Gather list of unique contigs
  contigs <- unique(manhat.base[,1])
  contigs <- contigs[which(!(is.na(contigs)))]

  #Create index data frame
  indexes <- as.data.frame(t(sapply(contigs,function(chr){
    return(c(chr,0,max(manhat.base[which(manhat.base[,1]==chr),2])))
  })))
  indexes$sum <- cumsum(indexes[,3])

  #Create new plotting dfs with modified coordinates & formatted values
  OR <- t(apply(OR,1,function(vals){
    #Quick fix to values
    vals <- as.vector(vals)
    vals[1] <- gsub("chr","",vals[1])
    vals[-1] <- as.numeric(vals[-1])

    #Get nearest plotting coordinate
    avgCoord <- mean(as.numeric(vals[2:3]))
    deltaCoords <- avgCoord-manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos
    nearestCoordIdx <- head(which(abs(deltaCoords)==min(abs(deltaCoords))),1)
    nearestCoord <- manhat.base[which(manhat.base$chr==as.character(vals[1])),]$pos[nearestCoordIdx]
    plottingCoord <- indexes[which(indexes[,1]==as.character(vals[1])),]$sum-indexes[which(indexes[,1]==as.character(vals[1])),3]+nearestCoord

    #Format ORs
    ORs <- as.numeric(as.vector(vals[-c(1:3)]))
    ORs[which(!is.infinite(ORs) & ORs>maxOR)] <- maxOR
    ORs[which(is.infinite(ORs) & ORs>0)] <- infOR
    ORs[which(is.na(ORs))] <- 0
    ORs <- log2(ORs)
    ORs[which(is.infinite(ORs) & ORs<0)] <- 0
    ORs[which(ORs<0)] <- 0

    #Return formatted values for plotting
    return(as.numeric(c(plottingCoord,ORs)))
  }))

  #Set height of contig boxes
  ymax <- log2(maxOR)

  #Prepare plotting window
  par(mar=c(0.5,2.5,0.2,0.2),bty="n")
  plot(x=c(0,4),y=c(0,1.02*ymax),type="n",
       xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Add gridlines
  abline(h=1:log2(maxOR),col=cols.CTRL[3])

  #Add axes
  y.at <- 0:ymax
  axis(2,at=y.at,labels=NA)
  axis(2,at=y.at,tick=F,line=-0.3,labels=2^y.at,las=2)
  axis(1,at=0:4,labels=NA,tck=0)

  #Plots points
  plot.df <- data.frame(c(2,4,14,24),c(23,13,23,36),
                        c(cols.GERM[1],cols.NEURO[1],
                          cols.SOMA[1],cols.CNCR[1]),
                        0.5:3.5)
  colnames(plot.df) <- c("start.idx","end.idx","col","pos")
  apply(plot.df,1,function(vals){
    #Clean up values
    vals <- as.vector(vals)
    index <- as.numeric(vals[1]):as.numeric(vals[2])

    #Get significant sites
    sig.idx <- apply(p[,index+2],1,function(p){
      if(any(p<=10^-8)){
        return(1)
      }else{
        return(0)
      }
    })
    sig.idx <- which(sig.idx>0)

    #Get mean ORs at significant sites
    means <- apply(OR[sig.idx,index],1,mean)

    #Plot swarm
    beeswarm(means,add=T,at=as.numeric(vals[4]),
             bg=vals[3],pch=21,
             method="swarm",corral="wrap",corralWidth=0.8)

    #Plot mean bar
    segments(x0=as.numeric(vals[4])-0.2,
             x1=as.numeric(vals[4])+0.2,
             y0=mean(means),y1=mean(means),
             lwd=3)
    segments(x0=as.numeric(vals[4])-0.2,
             x1=as.numeric(vals[4])+0.2,
             y0=mean(means),y1=mean(means),
             lwd=0.5,col="white")
  })
}


########################################
#####Generate manhattan plot for figure
########################################
#Plot mirrored manhattan
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegment_ORs.mirrorManhattan.pdf",sep=""),
    height=3.5,width=6.5)
mirrorManhattan(DEL.OR,DUP.OR,manhat.base)
dev.off()
#Plot mean del & dup ORs
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegment_meanORs.DEL.pdf",sep=""),
    height=2,width=2.75)
plotMedianORs(DEL.OR,DEL.p,manhat.base)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegment_meanORs.DUP.pdf",sep=""),
    height=2,width=2.75)
plotMedianORs(DUP.OR,DUP.p,manhat.base)
dev.off()


