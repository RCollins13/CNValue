#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Exploratory code for per-gene CNV burden scores

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
#Sample sizes (CTRL, GERM, CNCR)
nsamp <- c(38628,63629,10844)

#####Load required libraries
require(MASS)

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Read data
# GERM.E2.exonic <- sapply(c("CNV","DEL","DUP"),function(CNV){
#   dat <- read.table(paste(WRKDIR,"plot_data/perGene_burden/GERM_",CNV,
#                           "_E2_exonic.geneScore_data.txt",sep=""),
#                     header=T,comment.char="",sep="\t")
# })
GERM.N1.exonic <- sapply(c("CNV","DEL","DUP"),function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/perGene_burden/GERM_",CNV,
                          "_N1_exonic.geneScore_data.txt",sep=""),
                    header=T,comment.char="",sep="\t")
})
GERM.N1.wholegene <- sapply(c("CNV","DEL","DUP"),function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/perGene_burden/GERM_",CNV,
                          "_N1_wholegene.geneScore_data.txt",sep=""),
                    header=T,comment.char="",sep="\t")
})

#####Function to make scaled dotplot
plotDots <- function(df,CNV,nCTRL,nCASE,sim=10000,colP=T){
  #Gather plotting data
  dat <- data.frame(GERM.N1.exonic[,CNV]$control_CNV/nCTRL,
                    GERM.N1.exonic[,CNV]$case_CNV/nCASE)
  names(dat) <- c("CTRL","CASE")
  #Calculate Beta-Binomial ratio tests
  if(colP==T){
    # #CASE > CTRL
    # CASE.theta <- rbeta(sim,
    #                     (dat$CASE)+1,
    #                     (nCASE-dat$CASE)+1)
    # CTRL.theta <- rbeta(sim,
    #                     (dat$CTRL)+1,
    #                     (nCTRL-dat$CTRL)+1)
    # d.theta <- CASE.theta-CTRL.theta
    # thresholds <- quantile(d.theta,c(0.05/nrow(dat),0.05,0.95,1-(0.05/nrow(dat))))
    # #Test Beta-Binomial plot
    # plot(density(d.theta),
    #      xlab="CASE.theta - CTRL.theta",
    #      ylab="p(CASE.theta - CTRL.theta | nCASE, nCTRL)",
    #      main="Bayesian Posterior Probability of case CNV =/= control CNV",
    #      frame.plot=FALSE,cex.lab=1.5,lwd=3,yaxt="no")
    # abline(v=thresholds,lty=c(2,1,1,2),col=c("black","black","forestgreen","forestgreen"))
    # length(which(as.numeric(d.theta)>=as.numeric(thresholds[4])))
    # length(which(as.numeric(d.theta)>=as.numeric(thresholds[3])))
    # length(which(as.numeric(d.theta)<=as.numeric(thresholds[2])))
    # length(which(as.numeric(d.theta)<=as.numeric(thresholds[1])))
    #Fisher's exact apply
    pvals <- sapply(1:500,function(i){
      fisher.test(matrix(c(nCTRL-dat[i,1]*nCTRL,
                           nCASE-dat[i,2]*nCASE,
                           nCTRL,nCASE),
                         byrow=T,nrow=2))$p.value
    })
  }
  #Instantiate color vector
  colvect.sub <- sapply(pvals,function(p){
    if(p<=0.05){
      return("red")
    }else{
      return("black")
    }
  })
  colvect <- c(colvect.sub,rep(NA,nrow(dat)-length(sub)))
  #Set plotting parameters
  par(mar=c(4,4,0.5,0.5))
  #Determine axis maximum
  axMax <- max(dat)
  #Plot
  plot(dat$CTRL,dat$CASE,
       xlim=c(0,axMax),ylim=c(0,axMax),
       xlab="Control CNVs per Sample",
       ylab="Case CNVs per Sample",
       bg=colvect,pch=21)
  # #Overlay density plot
  # complete.dat <- dat[complete.cases(dat),]
  # contours <- kde2d(complete.dat[,1],complete.dat[,2],n=1000)
  # contour(contours,drawlabels=FALSE,nlevels=25,add=TRUE,
  #         col=adjustcolor(rev(rainbow(45)[1:32]),alpha=0.5))
  #X=Y line
  abline(0,1)
  #Linear regression line
  # abline(lm(nCASE ~ nCTRL,data=dat))
}

#####Calculate nominal and Bonferroni 95% CI bounds
# #Nominal
# GERM.nomP <- sapply(0:nsamp[1],function(nCTRL){
#   sapply(0:nsamp[2],function(nGERM){
#     fisher.test(matrix(c(nsamp[1]-nCTRL,nsamp[2]-nCASE,
#                          nCTRL,nCASE),
#                        byrow=T,nrow=2))$p.value
#   })
# })
# #Bonferroni
# GERM.bonP <- GERM.nomP*nrow(GERM[,1])

#####Test plots - case vs control
# #GERM E2 CNV
# plot(GERM.E2.exonic[,1]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,1]$case_CNV/nsamp[2])
# abline(0,1)
# #GERM E2 DEL
# plot(GERM.E2.exonic[,2]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,2]$case_CNV/nsamp[2])
# abline(0,1)
# #GERM E2 DUP
# plot(GERM.E2.exonic[,3]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,3]$case_CNV/nsamp[2])
# abline(0,1)
#GERM N1 exonic CNV
plotDots(GERM.N1.exonic,1,nsamp[1],nsamp[2],colP=F)
#GERM N1 exonic DEL
plotDots(GERM.N1.exonic,2,nsamp[1],nsamp[2],colP=F)
#GERM N1 exonic DUP
plotDots(GERM.N1.exonic,3,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene CNV
plotDots(GERM.N1.wholegene,1,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene DEL
plotDots(GERM.N1.wholegene,2,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene DUP
plotDots(GERM.N1.wholegene,3,nsamp[1],nsamp[2],colP=F)

#####Fisher test
dat <- head(as.data.frame(GERM.N1.exonic[,2]),500)
i <- which(dat$case_CNV==max(dat$case_CNV))
fisher.test(matrix(c(as.numeric(dat$control_CNV[i]),
                     as.numeric(dat$case_CNV[i]),
                     nsamp[1]-as.numeric(dat$control_CNV[i]),
                     nsamp[2]-as.numeric(dat$case_CNV[i])),
                   byrow=T,nrow=2))$p.value

