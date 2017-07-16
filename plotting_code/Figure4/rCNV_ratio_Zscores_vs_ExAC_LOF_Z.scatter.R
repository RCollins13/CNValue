#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate CNV Z-score rank vs ExAC scatterplots

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
    #Revise ExAC & RVIS percentiles based on subsetted gene list
    dat$lof_z_rank <- 100*rank(dat$lof_z_rank)/length(dat$lof_z_rank)
    dat$RVIS_pct <- 100*rank(dat$RVIS_pct)/length(dat$RVIS_pct)
    #Return data
    return(dat)
  })
}

######################################################
#####Helper function to plot ExAC vs rCNV scatterplots
######################################################
theta_ExAC_scatter <- function(df,yvar,color,
                               ylim=c(40,70),spar=0.35,
                               xaxis=T,yaxis=T){
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

  #Smooth response variable values
  yvar.stats.smoothed <- apply(yvar.stats,2,function(vals){
    smoothed <- smooth.spline(vals,spar=spar)
    return(data.frame(smoothed$x,smoothed$y))
  })

  #Prepare plot area
  par(mar=c(2,2,0.5,1))
  plot(x=c(0,101),y=ylim,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")

  #Add background shading & lines
  abline(h=seq(0,100,5),col=cols.CTRL[3])
  abline(h=50)
  rect(xleft=90,xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border=NA,col=adjustcolor("yellow",0.35))
  abline(v=90,lty=2)

  #Plot 95% CI
  polygon(x=c(yvar.stats.smoothed$V2[,1],rev(yvar.stats.smoothed$V3[,1])),
          y=c(yvar.stats.smoothed$V2[,2],rev(yvar.stats.smoothed$V3[,2])),
          border=NA,col=adjustcolor(color,0.4))

  #Plot observed percentiles
  points(x=1:100,y=yvar.stats[,1],lwd=2,pch=21,col=color,bg="white")

  #Plot smoothed trendline
  points(yvar.stats.smoothed$V1,type="l",lwd=4,col=color)

  #Add x-axis (if optioned)
  if(xaxis==T){
    axis(1,at=seq(0,100,5),labels=NA,tck=-0.01,col=cols.CTRL[1])
    axis(1,at=seq(0,100,10),labels=NA,tck=-0.02)
    axis(1,at=seq(0,100,20),tick=F,line=-0.8,cex.axis=1.3)
  }

  #Add y-axis (if optioned)
  if(yaxis==T){
    axis(2,at=seq(0,100,5),labels=NA,tck=-0.01,col=cols.CTRL[1])
    axis(2,at=seq(0,100,10),labels=NA,tck=-0.02)
    axis(2,at=seq(0,100,10),tick=F,line=-0.6,las=2,cex.axis=1.3)
  }
}

###############################################
#####Plot master theta pctile scatters for Fig4
###############################################
#Load data -- GERM/E4/exonic/minCNV=4
dat <- readData("GERM","E4","exonic",minCNV=1)

#Plot CNV
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_CNV.",
          "theta_RVIS_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[1]],yvar="RVIS",color="#333333",
                   ylim=c(40,65),xaxis=F)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_CNV.",
          "theta_constraint_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[1]],yvar="constraint",color="#333333",
                   ylim=c(40,65),xaxis=F,yaxis=F)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DEL.",
          "theta_RVIS_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[2]],yvar="RVIS",color="red",
                   ylim=c(40,65),xaxis=F)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DEL.",
          "theta_constraint_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[2]],yvar="constraint",color="red",
                   ylim=c(40,65),xaxis=F,yaxis=F)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DUP.",
          "theta_RVIS_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[3]],yvar="RVIS",color="blue",
                   ylim=c(40,65))
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DUP.",
          "theta_constraint_percentile_scatters.pdf",sep=""),
    width=4,height=2)
theta_ExAC_scatter(dat[[3]],yvar="constraint",color="blue",
                   ylim=c(40,65),yaxis=F)
dev.off()

