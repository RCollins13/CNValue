#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to make barplots of overlap between significant GERM & CNCR rCNV segments

####################################
#####Set parameters & load libraries
####################################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
require(beeswarm)
require(vioplot)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

##########################################
#####Helper function to read & format data
##########################################
readData <- function(pheno,CNV){
  dat <- read.table(paste(WRKDIR,"/plot_data/figure2/constrained_enrichments/",
                          pheno,".",CNV,".constrained_genes_count.txt",sep=""),header=F)
  colnames(dat) <- c("const","all")
  dat$frac <- dat$const/dat$all
  dat$frac[which(dat$all==0)] <- 0
  return(dat)
}

##########################################
#####Helper function to simulate null sets
##########################################
simNull <- function(dat,const=3264,all=19346,n=1000){
  #Simulate n sets
  sims <- matrix(unlist(lapply(1:n,function(i){
    #Set seed
    set.seed(i)
    #Hypergeometric simulation
    sim <- as.data.frame(t(sapply(dat$all,function(nGenes){
      if(nGenes>0){
        c(rhyper(1,const,all-const,nGenes),nGenes)
      }else{
        return(c(0,0))
      }
    })))
    colnames(sim) <- c("const","all")
    sim$frac <- sim$const/sim$all
    sim$frac[which(sim$all==0)] <- 0
    return(t(sim))
  })),ncol=3,byrow=T)

  #Clean up simulated data
  colnames(sims) <- c("const","all","frac")

  #Return simulated data
  return(as.data.frame(sims))
}

################################
#####Prep all data & simulations
################################
#Read data
ALL.DEL <- readData("ALL_PHENOS","DEL")
ALL.DUP <- readData("ALL_PHENOS","DUP")
GERM.DEL <- readData("GERM","DEL")
GERM.DUP <- readData("GERM","DUP")
CNCR.DEL <- readData("CNCR","DEL")
CNCR.DUP <- readData("CNCR","DUP")
#Simulate nulls
ALL.DEL.sim <- simNull(ALL.DEL)
ALL.DUP.sim <- simNull(ALL.DUP)
GERM.DEL.sim <- simNull(GERM.DEL)
GERM.DUP.sim <- simNull(GERM.DUP)
CNCR.DEL.sim <- simNull(CNCR.DEL)
CNCR.DUP.sim <- simNull(CNCR.DUP)

###########################################
#####Helper function: plot violins & swarms
###########################################
plotFracConst <- function(DEL.sim,DEL.obs,DUP.sim,DUP.obs){
  #Prepare plotting area
  par(mar=c(0.5,3,0.5,0.2),bty="n")
  plot(x=c(0,4),y=c(0,1),type="n",
       xaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Draw gridlines
  abline(h=seq(0,1,0.1),col=cols.CTRL[4])
  abline(h=seq(0,1,0.2),col=cols.CTRL[3])

  #Add axes
  axis(1,at=c(0,2,4),labels=NA)
  axis(2,at=seq(0,1,0.1),col=cols.CTRL[1],tck=-0.01,labels=NA)
  axis(2,at=seq(0,1,0.2),tck=-0.02,labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,las=2,line=-0.6,cex.axis=1.1,
       labels=paste(seq(0,100,20),"%",sep=""))

  #Plot violins of simulated nulls
  vioplot(DEL.sim$frac,at=0.5,col=cols.CTRL[1],add=T,pchMed=18,h=0.05)
  vioplot(DUP.sim$frac,at=2.5,col=cols.CTRL[1],add=T,pchMed=18,h=0.05)

  #Plot swarms of observed fractions
  beeswarm(DEL.obs$frac,add=T,at=1.5,bg="red",pch=21,cex=0.8,
           method="swarm",corral="wrap",corralWidth=0.8)
  beeswarm(DUP.obs$frac,add=T,at=3.5,bg="blue",pch=21,cex=0.8,
           method="swarm",corral="wrap",corralWidth=0.8)

  #Plot segments for means & medians
  segments(x0=c(1.2,3.2),x1=c(1.8,3.8),
           y0=mean(c(DEL.obs$frac,DUP.obs$frac)),
           y1=mean(c(DEL.obs$frac,DUP.obs$frac)),
           lwd=4)
  segments(x0=c(1.2,3.2),x1=c(1.8,3.8),
           y0=mean(c(DEL.obs$frac,DUP.obs$frac)),
           y1=mean(c(DEL.obs$frac,DUP.obs$frac)),
           lwd=1,col=c("red","blue"))
  segments(x0=c(1.2,3.2),x1=c(1.8,3.8),
           y0=mean(c(DEL.obs$frac,DUP.obs$frac)),
           y1=mean(c(DEL.obs$frac,DUP.obs$frac)),
           lwd=1,lty=2,col="white")

  # segments(x0=c(1.2,3.2),x1=c(1.8,3.8),
  #          y0=median(c(DEL.obs$frac,DUP.obs$frac)),
  #          y1=median(c(DEL.obs$frac,DUP.obs$frac)),
  #          lwd=2,lty=2)

  #M-W U-test of obs vs exp
  print(wilcox.test(DEL.obs$frac,DEL.sim$frac)$p.value)
  print(wilcox.test(DUP.obs$frac,DUP.sim$frac)$p.value)

  # #Add legend
  # legend("topright",legend=c("Mean","Median"),
  #        lty=c(1,3),lwd=2,bty="n")
}

#######################
#####Generate all plots
#######################
#Violin plots of constrained fractions
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_constrainedEnrichment.ALL.violins.pdf",sep=""),
    height=1.5,width=3)
plotFracConst(ALL.DEL.sim,ALL.DEL,ALL.DUP.sim,ALL.DUP)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_constrainedEnrichment.CNCR.violins.pdf",sep=""),
    height=1.5,width=3)
plotFracConst(CNCR.DEL.sim,CNCR.DEL,CNCR.DUP.sim,CNCR.DUP)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_constrainedEnrichment.GERM.violins.pdf",sep=""),
    height=1.5,width=3)
plotFracConst(GERM.DEL.sim,GERM.DEL,GERM.DUP.sim,GERM.DUP)
dev.off()





