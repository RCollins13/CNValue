#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to make four-way scatterplot of exonic/whole-gene del/dup Fisher's ORs per gene

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
require(plotrix)
require(MASS)

#################################################
#####Helper function to read data per combination
#################################################
readCNVs <- function(pheno,VF){
  #Load both DEL and DUP
  temp.dfs <- lapply(list("DEL","DUP"),function(CNV){
    #Load CNVs from both contexts
    temp.dfs <- lapply(list("exonic","wholegene"),function(context){
      #Load data & clean headers
      dat <- read.table(paste(WRKDIR,"plot_data/figure4/geneScore_results/",
                              pheno,"_",CNV,"_",VF,"_",context,
                              ".geneScore_stats.txt",sep=""),
                        comment.char="",header=T)
      colnames(dat)[1] <- "gene"

      #Retain three columns of interest
      dat.out <- data.frame(dat$gene,
                            dat$control_raw_CNV,
                            log2(dat$case_gt_control_Fisher_OR))

      #Only return genes with >0 control CNVs
      return(dat.out[which(dat.out$dat.control_raw_CNV>0),c(1,3)])
    })

    #Merge both dataframes
    dat.out <- merge(temp.dfs[[1]],temp.dfs[[2]],
                     by=1,suffixes=c("exonic","wholegene"))
    return(dat.out)
  })

  #Merge both dataframes
  dat.out <- merge(temp.dfs[[1]],temp.dfs[[2]],
                   by=1,suffixes=c("DEL","DUP"))
  return(dat.out)
}

######################################################
#####Helper function to plot scatter given two vectors
######################################################
plotScatter <- function(x,y,lim=6,plotlims=c(-6,6),ptcex=0.5,
                        xaxis=T,yaxis=T){
  #Prepare plotting area
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(x=plotlims,y=plotlims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")

  #Plot gridlines
  abline(h=-lim:lim,v=-lim:lim,col=cols.CTRL[3])
  abline(h=0,v=0,lwd=1.5,col=cols.CTRL[1])
  abline(0,1,col=cols.CTRL[1],lty=2)

  #Plot x-axis (if optioned)
  if(xaxis==T){
    axis(1,at=-lim:lim,col=cols.CTRL[1],
         labels=NA,tck=-0.01)
    axis(1,at=seq(-lim,lim,2),labels=NA,tck=-0.02)
    axis(1,at=seq(-lim,lim,2),tick=F,
         labels=c(paste("1/",2^seq(lim,1,-2),sep=""),
                  1,2^seq(2,lim,2)),line=-0.6,las=2,cex.axis=1.3)
  }

  #Plot y-axis (if optioned)
  if(yaxis==T){
    axis(2,at=-lim:lim,col=cols.CTRL[1],
         labels=NA,tck=-0.01)
    axis(2,at=seq(-lim,lim,2),labels=NA,tck=-0.02)
    axis(2,at=seq(-lim,lim,2),tick=F,
         labels=c(paste("1/",2^seq(lim,1,-2),sep=""),
                  1,2^seq(2,lim,2)),line=-0.6,las=2,cex.axis=1.3)
  }

  #Plot points
  points(x,y,pch=19,col=adjustcolor("black",alpha=0.2),cex=ptcex)

  #Plot trendline
  means.df <- data.frame(x,y)
  names(means.df) <- c("x","y")
  means.df <- means.df[which(!is.infinite(x) & ! is.infinite(y)),]
  means.df <- means.df[complete.cases(means.df),]
  abline(lm(y ~ x,data=means.df),lwd=2)

  #Plot density contours
  contours <- kde2d(means.df[,1],means.df[,2],n=200)
  contour(contours,drawlabels=FALSE,nlevels=30,add=TRUE,
          col=adjustcolor(rev(rainbow(50)[1:40]),alpha=0.5))

  #Compute stats (add manually in illustrator)
  print(cor(x,y,use="complete.obs",method="spearman"))
  print(lm(y ~ x,data=means.df)$coefficients[2])
  print(c(mean(2^x[which(!is.infinite((x)))],na.rm=T),
          mean(2^y[which(!is.infinite((y)))],na.rm=T)))
  print(mean(2^x[which(!is.infinite((x)))],na.rm=T)/mean(2^y[which(!is.infinite((y)))],na.rm=T))
  options(scipen=-1000)
  print(cor.test(x,y,method="spearman")$p.value)
  options(scipen=1000)
}

#Generate plots
# df <- readCNVs("GERM","E4")
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/FisherORs_exDEL_wgDEL.png",sep=""),
    width=3,height=3,units="in",res=1000)
plotScatter(df$log2.dat.case_gt_control_Fisher_OR.exonicDEL,
     df$log2.dat.case_gt_control_Fisher_OR.wholegeneDEL)
dev.off()
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/FisherORs_exDUP_wgDUP.png",sep=""),
    width=3,height=3,units="in",res=1000)
plotScatter(df$log2.dat.case_gt_control_Fisher_OR.exonicDUP,
     df$log2.dat.case_gt_control_Fisher_OR.wholegeneDUP,yaxis=F)
dev.off()
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/FisherORs_exDEL_exDUP.png",sep=""),
    width=3,height=3,units="in",res=1000)
plotScatter(df$log2.dat.case_gt_control_Fisher_OR.exonicDEL,
     df$log2.dat.case_gt_control_Fisher_OR.exonicDUP)
dev.off()
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/FisherORs_wgDEL_wgDUP.png",sep=""),
    width=3,height=3,units="in",res=1000)
plotScatter(df$log2.dat.case_gt_control_Fisher_OR.wholegeneDEL,
     df$log2.dat.case_gt_control_Fisher_OR.wholegeneDUP,yaxis=F)
dev.off()




