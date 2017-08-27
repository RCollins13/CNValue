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
  ymax <- log2(infOR)
  boxht <- 0.06*ymax

  #Prepare plotting window
  par(mar=c(0.2,2,0.2,0.2),bty="n")
  plot(x=c(-0.01*max(indexes[,4]),1.01*max(indexes[,4])),y=c(-1.1*(ymax+boxht),1.1*(ymax+boxht)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Add background shading & gridlines
  rect(xleft=indexes[,4]-indexes[,3],xright=indexes[,4],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=cols.CTRL[4:3],border=NA)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=-boxht,ytop=boxht,col="white",border=NA)
  abline(h=c((1:log2(maxOR))+boxht,(-1:-log2(maxOR))-boxht),col=cols.CTRL[2])
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=c(par("usr")[3],(log2(infOR)+boxht)+0.5*boxht),
       ytop=c((-log2(infOR)-boxht)-0.5*boxht,par("usr")[4]),
       border=NA,col="white")
  abline(h=c(ymax+boxht-0.5*boxht,ymax+boxht+0.5*boxht,-ymax-boxht-0.5*boxht,-ymax-boxht+0.5*boxht))
  rect(xleft=c(par("usr")[1],max(indexes[,4])),xright=c(0,par("usr")[2]),
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col="white",border=NA)

  #Add y-axis
  y.at <- 1:log2(maxOR)
  #Top y-axis
  axis(2,at=y.at+boxht,labels=NA)
  axis(2,at=y.at+boxht,tick=F,line=-0.3,labels=2^y.at,las=2)
  axis(2,at=ymax+boxht,labels=NA)
  axis(2,at=ymax+boxht,tick=F,line=-0.3,labels="Inf",las=2)
  #Bottom y-axis
  axis(2,at=-(y.at+boxht),labels=NA)
  axis(2,at=-(y.at+boxht),tick=F,line=-0.3,labels=2^y.at,las=2)
  axis(2,at=-ymax-boxht,labels=NA)
  axis(2,at=-ymax-boxht,tick=F,line=-0.3,labels="Inf",las=2)

  #Plots points
  apply(DEL.plot,1,function(vals){
    points(x=rep(vals[1],times=length(vals)-1),y=vals[-1]+boxht,
           pch=21,bg=cols.allPhenos[phenos.reorder],cex=0.5)
  })
  apply(DUP.plot,1,function(vals){
    points(x=rep(vals[1],times=length(vals)-1),y=-vals[-1]-boxht,
           pch=21,bg=cols.allPhenos[phenos.reorder],cex=0.5)
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

  #Add ticks where there's a significant
}

########################################
#####Generate manhattan plot for figure
########################################
#Plot
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegment_ORs.mirrorManhattan.png",sep=""),
    height=1200,width=1500,res=300)
mirrorManhattan(DEL.OR,DUP.OR,manhat.base)
dev.off()




