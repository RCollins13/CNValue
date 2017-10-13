#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to plot various tracks for SMARCA2 figure

####################################
#####Set parameters & load libraries
####################################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
nNDD <- 35693
nCTRL <- 38628
window.start <- 1250000
window.end <- 3250000
require(zoo)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#################################
#####Helper function to read data
#################################
readCNVs <- function(pheno,CNV,VF,filt,fileType="CNV"){
  #Read file
  if(fileType=="CNV"){
    in.path <- paste(WRKDIR,"/plot_data/SMARCA2/",pheno,".",CNV,".",
                     VF,".GRCh37.",filt,".bed.gz",sep="")
  }else{
    in.path <- paste(WRKDIR,"/plot_data/SMARCA2/",pheno,".",CNV,".",
                     VF,".GRCh37.",filt,".100bp_pileups.bed.gz",sep="")
  }
  dat <- read.table(in.path,header=F)

  #Add column names
  if(fileType=="CNV"){
    colnames(dat) <- c("chr","start","end","CNV_ID","CNV_type","phenos","PMID")
  }else{
    colnames(dat) <- c("chr","start","end","CNVs")
  }

  #Return data
  return(dat)
}

###################################################
#####Read 100bp pileups & normalize by sample sizes
###################################################
#Cases, E2, all
NDD.DEL.E2.all.pileup <- readCNVs("NDD","DEL","E2","all","pileup")
NDD.DEL.E2.all.pileup$CNVs.norm <- NDD.DEL.E2.all.pileup$CNVs/nNDD
NDD.DUP.E2.all.pileup <- readCNVs("NDD","DUP","E2","all","pileup")
NDD.DUP.E2.all.pileup$CNVs.norm <- NDD.DUP.E2.all.pileup$CNVs/nNDD
#Cases, E2, coding
NDD.DEL.E2.coding.pileup <- readCNVs("NDD","DEL","E2","coding","pileup")
NDD.DEL.E2.coding.pileup$CNVs.norm <- NDD.DEL.E2.coding.pileup$CNVs/nNDD
NDD.DUP.E2.coding.pileup <- readCNVs("NDD","DUP","E2","coding","pileup")
NDD.DUP.E2.coding.pileup$CNVs.norm <- NDD.DUP.E2.coding.pileup$CNVs/nNDD
#Cases, E2, haplosufficient
NDD.DEL.E2.haplosufficient.pileup <- readCNVs("NDD","DEL","E2","haplosufficient","pileup")
NDD.DEL.E2.haplosufficient.pileup$CNVs.norm <- NDD.DEL.E2.haplosufficient.pileup$CNVs/nNDD
NDD.DUP.E2.haplosufficient.pileup <- readCNVs("NDD","DUP","E2","haplosufficient","pileup")
NDD.DUP.E2.haplosufficient.pileup$CNVs.norm <- NDD.DUP.E2.haplosufficient.pileup$CNVs/nNDD
#Controls, E2, all
CTRL.DEL.E2.all.pileup <- readCNVs("CTRL","DEL","E2","all","pileup")
CTRL.DEL.E2.all.pileup$CNVs.norm <- CTRL.DEL.E2.all.pileup$CNVs/nCTRL
CTRL.DUP.E2.all.pileup <- readCNVs("CTRL","DUP","E2","all","pileup")
CTRL.DUP.E2.all.pileup$CNVs.norm <- CTRL.DUP.E2.all.pileup$CNVs/nCTRL
#Controls, E2, coding
CTRL.DEL.E2.coding.pileup <- readCNVs("CTRL","DEL","E2","coding","pileup")
CTRL.DEL.E2.coding.pileup$CNVs.norm <- CTRL.DEL.E2.coding.pileup$CNVs/nCTRL
CTRL.DUP.E2.coding.pileup <- readCNVs("CTRL","DUP","E2","coding","pileup")
CTRL.DUP.E2.coding.pileup$CNVs.norm <- CTRL.DUP.E2.coding.pileup$CNVs/nCTRL
#Cases, E2, haplosufficient
CTRL.DEL.E2.haplosufficient.pileup <- readCNVs("CTRL","DEL","E2","haplosufficient","pileup")
CTRL.DEL.E2.haplosufficient.pileup$CNVs.norm <- CTRL.DEL.E2.haplosufficient.pileup$CNVs/nCTRL
CTRL.DUP.E2.haplosufficient.pileup <- readCNVs("CTRL","DUP","E2","haplosufficient","pileup")
CTRL.DUP.E2.haplosufficient.pileup$CNVs.norm <- CTRL.DUP.E2.haplosufficient.pileup$CNVs/nCTRL

############################################################################
#####Function to plot mirrored pileups of case vs control, del & dup stacked
############################################################################
mirroredPileups <- function(case.del,case.dup,control.del,control.dup,
                            axis.scale=2/10000){
  #Prep dups for stacked plotting
  case.dup$stacked <- case.del$CNVs.norm+case.dup$CNVs.norm
  control.dup$stacked <- control.del$CNVs.norm+control.dup$CNVs.norm

  #Prep plotting area
  stackMax <- max(case.dup$stacked,control.dup$stacked)
  par(mar=c(0.4,1.8,0.4,0.2),bty="n")
  plot(x=c(window.start,window.end),y=1.05*c(-stackMax,stackMax),type="n",
       xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
  abline(h=0)

  #Add y-axis & gridlines
  axis(2,at=seq(0,stackMax,axis.scale),labels=NA)
  axis(2,at=seq(0,stackMax,axis.scale),tick=F,line=-0.2,las=2,
       labels=10000*seq(0,stackMax,axis.scale))
  abline(h=seq(0,stackMax,axis.scale/2),col=cols.CTRL[2])
  abline(h=seq(0,stackMax,axis.scale))
  axis(2,at=seq(0,-stackMax,-axis.scale),labels=NA)
  axis(2,at=seq(0,-stackMax,-axis.scale),tick=F,line=-0.2,las=2,
       labels=10000*seq(0,stackMax,axis.scale))
  abline(h=seq(0,-stackMax,-axis.scale/2),col=cols.CTRL[2])
  abline(h=seq(0,-stackMax,-axis.scale))

  #Plot stacked cases (positive y-axis)
  polygon(x=c(case.del$start,rev(case.del$start)),
          y=c(case.del$CNVs.norm,rep(0,nrow(case.del))),
          border=NA,col="red")
  polygon(x=c(case.del$start,rev(case.del$start)),
          y=c(case.dup$stacked,rev(case.del$CNVs.norm)),
          border=NA,col="blue")
  # points(x=case.dup$start,y=case.dup$stacked,type="l")

  #Plot stacked controls (negative y-axis)
  polygon(x=c(control.del$start,rev(control.del$start)),
          y=-c(control.del$CNVs.norm,rep(0,nrow(control.del))),
          border=NA,col="red")
  polygon(x=c(control.del$start,rev(control.del$start)),
          y=-c(control.dup$stacked,rev(control.del$CNVs.norm)),
          border=NA,col="blue")
  # points(x=control.dup$start,y=-control.dup$stacked,type="l")

  #Add midline
  abline(h=0,lwd=2)
}

#########################################
#####Function to plot odds ratios per bin
#########################################
plotORs <- function(case.del,case.dup,control.del,control.dup,q=0.95){
  #Collect sum of del+dup
  case.dup$stacked <- case.del$CNVs+case.dup$CNVs
  control.dup$stacked <- control.del$CNVs+control.dup$CNVs

  #Generate odds ratios
  del.ORs <- (case.del$CNVs/control.del$CNVs)/((nNDD-case.del$CNVs)/(nCTRL-control.del$CNVs))
  dup.ORs <- (case.dup$CNVs/control.dup$CNVs)/((nNDD-case.dup$CNVs)/(nCTRL-control.dup$CNVs))
  cnv.ORs <- (case.dup$stacked/control.dup$stacked)/((nNDD-case.dup$stacked)/(nCTRL-control.dup$stacked))

  #Round infinite odds ratios to max value
  allORs <- c(del.ORs,dup.ORs,cnv.ORs)
  maxVal <- floor(quantile(allORs[which(!is.infinite(allORs))],probs=q,na.rm=T))
  del.ORs[which(is.infinite(del.ORs) | del.ORs>maxVal)] <- maxVal
  dup.ORs[which(is.infinite(dup.ORs) | dup.ORs>maxVal)] <- maxVal
  cnv.ORs[which(is.infinite(cnv.ORs) | cnv.ORs>maxVal)] <- maxVal
  del.ORs[which(is.na(del.ORs))] <- 1
  dup.ORs[which(is.na(dup.ORs))] <- 1
  cnv.ORs[which(is.na(cnv.ORs))] <- 1

  #Smooth odds ratios
  del.ORs.smooth <- rollmean(del.ORs,k=100)
  dup.ORs.smooth <- rollmean(dup.ORs,k=100)
  cnv.ORs.smooth <- rollmean(cnv.ORs,k=100)

  #Prep plotting area
  par(mar=c(0.4,1.8,0.4,0.2),bty="n")
  plot(x=c(window.start,window.end),y=c(1,maxVal),type="n",
       xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
  abline(h=c(1,maxVal))

  #Plot smoothed odds ratios
  points(x=case.del$start[-c(1:49,(length(del.ORs)-49):length(del.ORs))],y=del.ORs.smooth,
         type="l",lwd=3,col="red")
  points(x=case.dup$start[-c(1:49,(length(dup.ORs)-49):length(dup.ORs))],y=dup.ORs.smooth,
         type="l",lwd=3,col="blue")
  # points(x=case.del$start[-c(1:49,(length(cnv.ORs)-49):length(cnv.ORs))],y=cnv.ORs.smooth,
  #        type="l",lwd=3,col="black")

  #Add axis
  axis(2,at=c(1,floor(maxVal)),labels=NA)
  axis(2,at=c(1,floor(maxVal)),tick=F,
       line=-0.2,las=2,labels=c(1,floor(maxVal)))
}

################
#####Plot tracks
################
#Prepare plotting area -- all CNVs
pdf(paste(WRKDIR,"/rCNV_map_paper/Figures/SMARCA2/SMARCA2_CNV_tracks.all_CNVs.pdf",sep=""),
    height=2,width=12)
layout(matrix(c(1,2),nrow=2),heights=c(5,2))
#Mirrored pileup -- all CNVs
mirroredPileups(NDD.DEL.E2.all.pileup,
                NDD.DUP.E2.all.pileup,
                CTRL.DEL.E2.all.pileup,
                CTRL.DUP.E2.all.pileup,
                axis.scale=8/10000)
#Plot ORs -- all CNVs
plotORs(NDD.DEL.E2.all.pileup,
        NDD.DUP.E2.all.pileup,
        CTRL.DEL.E2.all.pileup,
        CTRL.DUP.E2.all.pileup,
        q=0.90)
#Finish plotting
dev.off()

#Prepare plotting area -- coding CNVs
pdf(paste(WRKDIR,"/rCNV_map_paper/Figures/SMARCA2/SMARCA2_CNV_tracks.coding_CNVs.pdf",sep=""),
    height=2,width=12)
layout(matrix(c(1,2),nrow=2),heights=c(5,2))
#Mirrored pileup -- coding CNVs
mirroredPileups(NDD.DEL.E2.coding.pileup,
                NDD.DUP.E2.coding.pileup,
                CTRL.DEL.E2.coding.pileup,
                CTRL.DUP.E2.coding.pileup,
                axis.scale=8/10000)
#Plot ORs -- coding CNVs
plotORs(NDD.DEL.E2.coding.pileup,
        NDD.DUP.E2.coding.pileup,
        CTRL.DEL.E2.coding.pileup,
        CTRL.DUP.E2.coding.pileup,
        q=0.90)
#Finish plotting
dev.off()

#Prepare plotting area -- haplosufficient CNVs
pdf(paste(WRKDIR,"/rCNV_map_paper/Figures/SMARCA2/SMARCA2_CNV_tracks.haplosufficient_CNVs.pdf",sep=""),
    height=2,width=12)
layout(matrix(c(1,2),nrow=2),heights=c(5,2))
#Mirrored pileup -- haplosufficient CNVs
mirroredPileups(NDD.DEL.E2.haplosufficient.pileup,
                NDD.DUP.E2.haplosufficient.pileup,
                CTRL.DEL.E2.haplosufficient.pileup,
                CTRL.DUP.E2.haplosufficient.pileup,
                axis.scale=4/10000)
#Plot ORs -- haplosufficient CNVs
plotORs(NDD.DEL.E2.haplosufficient.pileup,
        NDD.DUP.E2.haplosufficient.pileup,
        CTRL.DEL.E2.haplosufficient.pileup,
        CTRL.DUP.E2.haplosufficient.pileup,
        q=0.99)
dev.off()


