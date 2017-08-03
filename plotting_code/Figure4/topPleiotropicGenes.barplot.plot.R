#!/usr/bin/env R

#rCNV Map Project
#Summer 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate herringbone barplot of most recurrently mutated (pleiotropic) genes

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

##########################################
#####Helper function to read & format data
##########################################
readDat <- function(path){
  #Read data
  dat <- read.table(path,header=F)

  #Convert group counts to pct
  dat <- cbind(dat,
               t(apply(dat[,2:3],1,function(vals){return(vals/22)})),
               t(apply(dat[,4:5],1,function(vals){return(vals/13)})))
  dat$max <- apply(dat[,6:9],1,max)
  dat$max.GERM <- apply(dat[,6:7],1,max)
  dat$max.CNCR <- apply(dat[,8:9],1,max)
  dat$avg <- apply(dat[,6:9],1,mean)

  #Reformat column names
  names(dat) <- c("gene","GERM.DEL","GERM.DUP","CNCR.DEL","CNCR.DUP",
                  "GERM.DEL.pct","GERM.DUP.pct","CNCR.DEL.pct","CNCR.DUP.pct",
                  "max","max.GERM","max.CNCR","mean")
  return(dat)
}

##############
#####Read data
##############
dat <- readDat(paste(WRKDIR,"/plot_data/figure4/sigGenes_min6_countsPerGene_byGroup.txt",sep=""))
dat.orderGERM <- dat[order(-dat$max.GERM,-dat$max.CNCR,dat$gene),]
dat.orderCNCR <- dat[order(-dat$max.CNCR,-dat$max.GERM,dat$gene),]

#####Helper funnction to plot horizontal herringbone plots
herringbone <- function(df,N=25,cMar=0.3,rev=F,
                        ytop=100,ybottom=100,
                        highlight=NULL){
  #Invert df & set margins if optioned
  if(rev==T){
    df <- head(df,N)[N:1,]
    mar <- c(0.5,0.5,0.5,2.5)
    axSide <- 4
  }else{
    mar <- c(0.5,2.5,0.5,0.5)
    axSide <- 2
  }

  #Prepare highlights
  fills <- rep(cols.CTRL[3],N)
  fonts <- rep(3,N)
  borders <- rep(NA,N)
  if(!is.null(highlight)){
    fills[highlight] <- "yellow"
    fonts[highlight] <- 4
    borders[highlight] <- "black"
  }

  #Prepare plot area
  par(mar=mar,bty="n")
  plot(x=c(0,N),y=c(-1-cMar,1+cMar),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")

  #Add plot background
  rect(xleft=(0:(N-1))+0.1,xright=(1:N)-0.1,
       ybottom=-cMar,ytop=cMar,
       border=borders,col=fills)
  abline(h=seq(cMar,cMar+(ytop/100),0.25),col=cols.CTRL[4])
  abline(h=-seq(cMar,cMar+(ybottom/100),0.25),col=cols.CTRL[4])

  #Iterate and plot bars
  sapply(1:N,function(i){
    rect(xleft=c(i-0.8,i-0.5,i-0.8,i-0.5),
         xright=c(i-0.5,i-0.2,i-0.5,i-0.2),
         ybottom=c(cMar,cMar,-df[i,8:9]-cMar),
         ytop=c(df[i,6:7]+cMar,-cMar,-cMar),
         col=c("red","blue","red","blue"),lwd=0.5)
  })

  #Add cleanup lines
  abline(h=c(-cMar,cMar))

  #Add gene symbols
  sapply(1:N,function(i){
    text(x=i-0.5,y=0,srt=90,font=fonts[i],labels=df[i,1])
  })

  #Add top y-axis
  axis(axSide,at=seq(cMar,(ytop/100)+cMar,0.125),
       labels=NA,tck=-0.015,col=cols.GERM[3])
  axis(axSide,at=seq(cMar,(ytop/100)+cMar,0.25),
       labels=NA,tck=-0.025,col=cols.GERM[1])
  axis(axSide,at=seq(cMar,(ytop/100)+cMar,0.25),tick=F,
       labels=paste(seq(0,ytop,25),"%",sep=""),
       las=2,line=-0.6,cex.axis=0.9)

  #Add bottom y-axis
  axis(axSide,at=-seq(cMar,(ybottom/100)+cMar,0.125),
       labels=NA,tck=-0.015,col=cols.CNCR[3])
  axis(axSide,at=-seq(cMar,(ybottom/100)+cMar,0.25),
       labels=NA,tck=-0.025,col=cols.CNCR[1])
  axis(axSide,at=-seq(cMar,(ybottom/100)+cMar,0.25),tick=F,
       labels=paste(seq(0,ybottom,25),"%",sep=""),
       las=2,line=-0.6,cex.axis=0.9)
}

###################
#####Generate plots
###################
#GERM
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/topPleiotropicGenes_GERM.pdf",sep=""),
    width=9,height=3)
herringbone(dat.orderGERM,N=40,cMar=0.5,ybottom=50,
            highlight=c(6,12,23,24,31,40))
dev.off()
#CNCR
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/topPleiotropicGenes_CNCR.pdf",sep=""),
    width=9,height=3)
herringbone(dat.orderCNCR,N=40,cMar=0.5,ytop=50,
            highlight=c(23,39))
dev.off()





