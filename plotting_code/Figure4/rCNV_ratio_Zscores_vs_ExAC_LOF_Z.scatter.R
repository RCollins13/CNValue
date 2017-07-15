#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate CNV Z-score rank vs ExAC  scatterplots

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

##################################################################################
#####Helper function to read data, & match with ExAC scores (based on gene symbol)
##################################################################################
#Read ExAC scores
ExAC.constraint <- read.table(paste(WRKDIR,"plot_data/figure4/ExAC_LoF_constraint.txt",sep=""),
                     header=F)
colnames(ExAC.constraint) <- c("gene","ObsExp","lof_z","pLI")
ExAC.constraint$lof_z_rank <- 100*rank(ExAC.constraint$lof_z)/length(ExAC.constraint$lof_z)
ExAC.RVIS <- read.table(paste(WRKDIR,"plot_data/figure4/ExAC_RVIS.txt",sep=""),
                        header=F)
colnames(ExAC.RVIS) <- c("gene","RVIS","RVIS_pct")
#Read CNV data & match w/ExAC Z-scores
readData <- function(pheno,VF,context,minCNV=2){
  lapply(list("CNV","DEL","DUP"),function(CNV){
    #Read table
    dat <- read.table(paste(WRKDIR,"plot_data/figure4/geneScore_results/",
                            pheno,"_",CNV,"_",VF,"_",context,
                            ".geneScore_stats.txt",sep=""),
                      comment.char="",header=T)
    #Exclude genes with fewer than minCNV total CNVs
    dat <- dat[which(dat$total_raw_CNV>=minCNV),]
    #Replace first column header
    names(dat)[1] <- "gene"
    #Match w/ExAC
    dat <- merge(dat,ExAC.constraint)
    #Match w/RVIS
    dat <- merge(dat,ExAC.RVIS)
    #Return data
    return(dat)
  })
}

######################################################
#####Helper function to plot ExAC vs rCNV scatterplots
######################################################
theta_ExAC_scatter <- function(df,yvar,color){
  #Compute rCNV z-score percentiles
  rCNV.cents <- quantile(x=df$theta_Zscore,probs=seq(0,1,0.01))

  #Compute response variable mean & 95% CIs per centile
  yvar.stats <- as.data.frame(t(sapply(1:100,function(q){
    #Extract relevant values
    if(yvar=="constraint"){
      vals <- df[which(df$theta_Zscore>=rCNV.cents[q]
                       & df$theta_Zscore<rCNV.cents[q+1]),]$lof_z_rank
    }else{
      vals <- 100-df[which(df$theta_Zscore>=rCNV.cents[q]
                       & df$theta_Zscore<rCNV.cents[q+1]),]$RVIS_pct
    }

    #Compute stats & return
    CI.margin <- 1.96*std.error(vals)
    stats <- c(mean(vals),
               mean(vals)-CI.margin,
               mean(vals)+CI.margin)
    return(stats)
  })))
}


#Test percentile plots - theta Z vs ExAC lof Z
dat <- readData("GERM","E4","exonic",minCNV=3)


#DEL

DEL.means <- sapply(1:100,function(q){
  mean(c(DEL.cents[q],DEL.cents[q+1]))
})
DEL.lof_z <- sapply(1:100,function(q){
  df <- dat[[2]]

})
plot(DEL.means,DEL.lof_z,xlim=c(-2,2))
DEL.lof_z_pct <- sapply(1:100,function(q){
  df <- dat[[2]]
  mean(df[which(df$theta_Zscore>=DEL.cents[q] & df$theta_Zscore<DEL.cents[q+1]),]$lof_z_rank)
})
plot(DEL.means,DEL.RVIS_pct,xlim=c(-2,2))
###############
plot(DEL.lof_z_pct,lwd=2,pch=21,col="red",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
abline(h=50)
points(smooth.spline(1:100,DEL.lof_z_pct,spar=0.4),lwd=3,col="red",type="l")
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############
DEL.RVIS <- sapply(1:100,function(q){
  df <- dat[[2]]
  mean(df[which(df$theta_Zscore>=DEL.cents[q] & df$theta_Zscore<DEL.cents[q+1]),]$RVIS)
})
plot(DEL.means,DEL.RVIS,xlim=c(-2,2))
DEL.RVIS_pct <- sapply(1:100,function(q){
  df <- dat[[2]]
  mean(df[which(df$theta_Zscore>=DEL.cents[q] & df$theta_Zscore<DEL.cents[q+1]),]$RVIS_pct)
})
plot(DEL.means,DEL.RVIS_pct,xlim=c(-2,2))
###############
plot(100-DEL.RVIS_pct,lwd=2,pch=21,col="red",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
abline(h=50)
points(smooth.spline(1:100,100-DEL.RVIS_pct,spar=0.4),lwd=3,col="red",type="l")
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############


#DUP
DUP.cents <- quantile(x=dat[[3]]$theta_Zscore,probs=seq(0,1,0.01))
DUP.means <- sapply(1:100,function(q){
  mean(c(DUP.cents[q],DUP.cents[q+1]))
})
DUP.lof_z <- sapply(1:100,function(q){
  df <- dat[[3]]
  mean(df[which(df$theta_Zscore>=DUP.cents[q] & df$theta_Zscore<DUP.cents[q+1]),]$lof_z)
})
plot(DUP.means,DUP.lof_z,xlim=c(-5,5))
DUP.lof_z_pct <- sapply(1:100,function(q){
  df <- dat[[3]]
  mean(df[which(df$theta_Zscore>=DUP.cents[q] & df$theta_Zscore<DUP.cents[q+1]),]$lof_z_rank)
})
plot(DUP.means,DUP.lof_z_pct,xlim=c(-2,2))
###############
plot(DUP.lof_z_pct,lwd=2,pch=21,col="blue",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
abline(h=50)
points(smooth.spline(1:100,DUP.lof_z_pct,spar=0.4),lwd=3,col="blue",type="l")
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############
DUP.RVIS <- sapply(1:100,function(q){
  df <- dat[[3]]
  mean(df[which(df$theta_Zscore>=DUP.cents[q] & df$theta_Zscore<DUP.cents[q+1]),]$RVIS)
})
plot(DUP.means,DUP.RVIS,xlim=c(-2,2))
DUP.RVIS_pct <- sapply(1:100,function(q){
  df <- dat[[3]]
  mean(df[which(df$theta_Zscore>=DUP.cents[q] & df$theta_Zscore<DUP.cents[q+1]),]$RVIS_pct)
})
plot(DUP.means,DUP.RVIS_pct,xlim=c(-2,2))
###############
plot(100-DUP.RVIS_pct,lwd=2,pch=21,col="blue",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
points(smooth.spline(1:100,100-DUP.RVIS_pct,spar=0.4),lwd=3,col="blue",type="l")
abline(h=50)
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############


#CNV
CNV.cents <- quantile(x=dat[[1]]$theta_Zscore,probs=seq(0,1,0.01))
CNV.means <- sapply(1:100,function(q){
  mean(c(CNV.cents[q],CNV.cents[q+1]))
})
CNV.lof_z <- sapply(1:100,function(q){
  df <- dat[[1]]
  mean(df[which(df$theta_Zscore>=CNV.cents[q] & df$theta_Zscore<CNV.cents[q+1]),]$lof_z)
})
plot(CNV.means,CNV.lof_z,xlim=c(-2,2))
CNV.lof_z_pct <- sapply(1:100,function(q){
  df <- dat[[1]]
  mean(df[which(df$theta_Zscore>=CNV.cents[q] & df$theta_Zscore<CNV.cents[q+1]),]$lof_z_rank)
})
plot(CNV.means,CNV.lof_z_pct,xlim=c(-2,2))
###############
plot(CNV.lof_z_pct,lwd=2,pch=21,col="black",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
abline(h=50)
points(smooth.spline(1:100,CNV.lof_z_pct,spar=0.4),lwd=3,col="black",type="l")
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############
CNV.RVIS <- sapply(1:100,function(q){
  df <- dat[[1]]
  mean(df[which(df$theta_Zscore>=CNV.cents[q] & df$theta_Zscore<CNV.cents[q+1]),]$RVIS)
})
plot(CNV.means,CNV.RVIS,xlim=c(-2,2))
CNV.RVIS_pct <- sapply(1:100,function(q){
  df <- dat[[1]]
  mean(df[which(df$theta_Zscore>=CNV.cents[q] & df$theta_Zscore<CNV.cents[q+1]),]$RVIS_pct)
})
plot(CNV.means,CNV.RVIS_pct,xlim=c(-2,2))
###############
plot(100-CNV.RVIS_pct,lwd=2,pch=21,col="black",
     xlab="rCNV Burden Percentile",ylim=c(45,65),
     panel.first=c(rect(xleft=90,xright=par("usr")[2],
                        ybottom=par("usr")[3],ytop=par("usr")[4],
                        col="yellow",border=NA)))
abline(h=50)
points(smooth.spline(1:100,100-CNV.RVIS_pct,spar=0.4),lwd=3,col="black",type="l")
abline(v=90,lty=2)
# abline(v=95,lty=3)
###############



