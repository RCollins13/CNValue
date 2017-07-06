#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate dotplots of rCNV burdens for gene set categories

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
phenos.GERM <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
                 "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
                 "MSC","EE","INT","EMI")
phenos.CNCR <- c("CNCR","CGEN","CSKN","CGST","CRNL",
                 "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(plotrix)

#####Read data
OR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                       "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
#Subset by GERM and CNCR
GERM <- OR[,2:23]
CNCR <- OR[,24:36]

#####Gather means & 95% CIs for ALL, GERM and CNCR
ALL.mean <- mean(log2(apply(OR[,-1],2,mean,na.rm=T)),na.rm=T)
ALL.CI <- 1.96*std.error(log2(apply(OR[,-1],2,mean,na.rm=T)),na.rm=T)
GERM.mean <- mean(log2(apply(GERM,2,mean,na.rm=T)),na.rm=T)
GERM.CI <- 1.96*std.error(log2(apply(GERM,2,mean,na.rm=T)),na.rm=T)
CNCR.mean <- mean(log2(apply(CNCR,2,mean,na.rm=T)),na.rm=T)
CNCR.CI <- 1.96*std.error(log2(apply(CNCR,2,mean,na.rm=T)),na.rm=T)

#####Simulate expected gene set if normally distributed around OR=1
EXP.mean <- mean(rnorm(mean=0,n=nrow(OR),
                       sd=sd(log2(apply(OR[,-1],2,mean,na.rm=T)),na.rm=T)),na.rm=T)
EXP.CI <- 1.96*std.error(rnorm(mean=0,n=nrow(OR),
                               sd=sd(log2(apply(OR[,-1],2,mean,na.rm=T)),na.rm=T)))

#####Gather tissue-agnostic means & CIs for ALL, GERM, and CNCR
#Haplosufficient genes
ALL.Haplosuff.avg.mean <- mean(log2(apply(OR[3,-1],2,mean,na.rm=T)),na.rm=T)
ALL.Haplosuff.avg.CI <- 1.96*std.error(log2(apply(OR[3,-1],2,mean,na.rm=T)),na.rm=T)
GERM.Haplosuff.avg.mean <- mean(log2(apply(GERM[3,],2,mean,na.rm=T)),na.rm=T)
GERM.Haplosuff.avg.CI <- 1.96*std.error(log2(apply(GERM[3,],2,mean,na.rm=T)),na.rm=T)
CNCR.Haplosuff.avg.mean <- mean(log2(apply(CNCR[3,],2,mean,na.rm=T)),na.rm=T)
CNCR.Haplosuff.avg.CI <- 1.96*std.error(log2(apply(CNCR[3,],2,mean,na.rm=T)),na.rm=T)
#Not haplosufficient genes
ALL.NotHaplosuff.avg.mean <- mean(log2(apply(OR[4,-1],2,mean,na.rm=T)),na.rm=T)
ALL.NotHaplosuff.avg.CI <- 1.96*std.error(log2(apply(OR[4,-1],2,mean,na.rm=T)),na.rm=T)
GERM.NotHaplosuff.avg.mean <- mean(log2(apply(GERM[4,],2,mean,na.rm=T)),na.rm=T)
GERM.NotHaplosuff.avg.CI <- 1.96*std.error(log2(apply(GERM[4,],2,mean,na.rm=T)),na.rm=T)
CNCR.NotHaplosuff.avg.mean <- mean(log2(apply(CNCR[4,],2,mean,na.rm=T)),na.rm=T)
CNCR.NotHaplosuff.avg.CI <- 1.96*std.error(log2(apply(CNCR[4,],2,mean,na.rm=T)),na.rm=T)
#Constrained genes
ALL.Const.avg.mean <- mean(log2(apply(OR[5,-1],2,mean,na.rm=T)),na.rm=T)
ALL.Const.avg.CI <- 1.96*std.error(log2(apply(OR[5,-1],2,mean,na.rm=T)),na.rm=T)
GERM.Const.avg.mean <- mean(log2(apply(GERM[5,],2,mean,na.rm=T)),na.rm=T)
GERM.Const.avg.CI <- 1.96*std.error(log2(apply(GERM[5,],2,mean,na.rm=T)),na.rm=T)
CNCR.Const.avg.mean <- mean(log2(apply(CNCR[5,],2,mean,na.rm=T)),na.rm=T)
CNCR.Const.avg.CI <- 1.96*std.error(log2(apply(CNCR[5,],2,mean,na.rm=T)),na.rm=T)
#Highly constrained genes
ALL.HiConst.avg.mean <- mean(log2(apply(OR[6,-1],2,mean,na.rm=T)),na.rm=T)
ALL.HiConst.avg.CI <- 1.96*std.error(log2(apply(OR[6,-1],2,mean,na.rm=T)),na.rm=T)
GERM.HiConst.avg.mean <- mean(log2(apply(GERM[6,],2,mean,na.rm=T)),na.rm=T)
GERM.HiConst.avg.CI <- 1.96*std.error(log2(apply(GERM[6,],2,mean,na.rm=T)),na.rm=T)
CNCR.HiConst.avg.mean <- mean(log2(apply(CNCR[6,],2,mean,na.rm=T)),na.rm=T)
CNCR.HiConst.avg.CI <- 1.96*std.error(log2(apply(CNCR[6,],2,mean,na.rm=T)),na.rm=T)
#Extremely constrained genes
ALL.ExtConst.avg.mean <- mean(log2(apply(OR[7,-1],2,mean,na.rm=T)),na.rm=T)
ALL.ExtConst.avg.CI <- 1.96*std.error(log2(apply(OR[7,-1],2,mean,na.rm=T)),na.rm=T)
GERM.ExtConst.avg.mean <- mean(log2(apply(GERM[7,],2,mean,na.rm=T)),na.rm=T)
GERM.ExtConst.avg.CI <- 1.96*std.error(log2(apply(GERM[7,],2,mean,na.rm=T)),na.rm=T)
CNCR.ExtConst.avg.mean <- mean(log2(apply(CNCR[7,],2,mean,na.rm=T)),na.rm=T)
CNCR.ExtConst.avg.CI <- 1.96*std.error(log2(apply(CNCR[7,],2,mean,na.rm=T)),na.rm=T)

#####Load pheno-geneset matching table
matches <- read.table(paste(WRKDIR,"rCNVmap/misc/phenotype_geneSet_matches.list",sep=""),header=F)

#####Get average of categories of interest across all tissues
#OMIM-HPO
ALL.OMIM.avg.mean <- mean(log2(apply(OR[28:61,-1],2,mean,na.rm=T)),na.rm=T)
ALL.OMIM.avg.CI <- 1.96*std.error(log2(apply(OR[28:61,-1],2,mean,na.rm=T)),na.rm=T)
GERM.OMIM.avg.mean <- mean(log2(apply(GERM[28:61,],2,mean,na.rm=T)),na.rm=T)
GERM.OMIM.avg.CI <- 1.96*std.error(log2(apply(GERM[28:61,],2,mean,na.rm=T)),na.rm=T)
CNCR.OMIM.avg.mean <- mean(log2(apply(CNCR[28:61,],2,mean,na.rm=T)),na.rm=T)
CNCR.OMIM.avg.CI <- 1.96*std.error(log2(apply(CNCR[28:61,],2,mean,na.rm=T)),na.rm=T)
#Highly expressed
HighExpr.idx <- c(67,74,81,88,95,102,109,116,123,130,137,143,150)
ALL.HighExpr.avg.mean <- mean(log2(apply(OR[HighExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
ALL.HighExpr.avg.CI <- 1.96*std.error(log2(apply(OR[HighExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
GERM.HighExpr.avg.mean <- mean(log2(apply(GERM[HighExpr.idx,],2,mean,na.rm=T)),na.rm=T)
GERM.HighExpr.avg.CI <- 1.96*std.error(log2(apply(GERM[HighExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.HighExpr.avg.mean <- mean(log2(apply(CNCR[HighExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.HighExpr.avg.CI <- 1.96*std.error(log2(apply(CNCR[HighExpr.idx,],2,mean,na.rm=T)),na.rm=T)
#Very highly expressed
VeryExpr.idx <- c(72,79,86,93,100,107,114,121,128,135,141,148,155)
ALL.VeryExpr.avg.mean <- mean(log2(apply(OR[VeryExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
ALL.VeryExpr.avg.CI <- 1.96*std.error(log2(apply(OR[VeryExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
GERM.VeryExpr.avg.mean <- mean(log2(apply(GERM[VeryExpr.idx,],2,mean,na.rm=T)),na.rm=T)
GERM.VeryExpr.avg.CI <- 1.96*std.error(log2(apply(GERM[VeryExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.VeryExpr.avg.mean <- mean(log2(apply(CNCR[VeryExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.VeryExpr.avg.CI <- 1.96*std.error(log2(apply(CNCR[VeryExpr.idx,],2,mean,na.rm=T)),na.rm=T)
#Specifically expressed
SpecExpr.idx <- c(71,78,85,92,99,106,113,120,127,134,147,154)
ALL.SpecExpr.avg.mean <- mean(log2(apply(OR[SpecExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
ALL.SpecExpr.avg.CI <- 1.96*std.error(log2(apply(OR[SpecExpr.idx,-1],2,mean,na.rm=T)),na.rm=T)
GERM.SpecExpr.avg.mean <- mean(log2(apply(GERM[SpecExpr.idx,],2,mean,na.rm=T)),na.rm=T)
GERM.SpecExpr.avg.CI <- 1.96*std.error(log2(apply(GERM[SpecExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.SpecExpr.avg.mean <- mean(log2(apply(CNCR[SpecExpr.idx,],2,mean,na.rm=T)),na.rm=T)
CNCR.SpecExpr.avg.CI <- 1.96*std.error(log2(apply(CNCR[SpecExpr.idx,],2,mean,na.rm=T)),na.rm=T)

#####Get value of categories matched between phenos & genos
matched <- lapply(list("OMIM","HighExpr","VeryExpr","SpecExpr"),function(category){
  #Get all matched ORs
  OR.vals <- sapply(unique(matches[which(matches[,2]==category),1]),function(pheno){
    vals <- sapply(unique(matches[which(matches[,2]==category & matches[,1]==pheno),3]),function(geneset){
      return(OR[which(OR[,1]==geneset),which(phenos==pheno)+1])
    })
    return(mean(log2(vals),na.rm=T))
  })

  #Calculate ALL, GERM, and CNCR means & 95% CIs
  ALL.m <- mean(OR.vals,na.rm=T)
  ALL.CI <- 1.96*std.error(OR.vals,na.rm=T)
  GERM.m <- mean(OR.vals[which(names(OR.vals) %in% phenos.GERM)],na.rm=T)
  GERM.CI <- 1.96*std.error(OR.vals[which(names(OR.vals) %in% phenos.GERM)],na.rm=T)
  CNCR.m <- mean(OR.vals[which(names(OR.vals) %in% phenos.CNCR)],na.rm=T)
  CNCR.CI <- 1.96*std.error(OR.vals[which(names(OR.vals) %in% phenos.CNCR)],na.rm=T)

  #Return all data
  return(c(ALL.m,ALL.CI,GERM.m,GERM.CI,CNCR.m,CNCR.CI))
})

#####Get values for constraint-expression combination categories
combs.est <- read.table(paste(WRKDIR,"plot_data/figure3/constraint_expression_combinations/",
                       "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
combs.est[,-1] <- log2(combs.est[,-1])
combs.ALL.m <- apply(combs.est[,-1],1,mean,na.rm=T)
combs.ALL.CI <- 1.96*apply(combs.est[,-1],1,std.error,na.rm=T)
combs.GERM.m <- apply(combs.est[,2:23],1,mean,na.rm=T)
combs.GERM.CI <- 1.96*apply(combs.est[,2:23],1,std.error,na.rm=T)
combs.CNCR.m <- apply(combs.est[,24:36],1,mean,na.rm=T)
combs.CNCR.CI <- 1.96*apply(combs.est[,24:36],1,std.error,na.rm=T)

#####Prepare plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/geneSet_categories_ORs.dotplot.pdf",sep=""),
    width=2.5,height=4.5)
par(mar=c(2.5,1,0.5,0.75),bty="n")
plot(x=log2(c(1,2.5)),y=c(0,-21),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")

#####Draw gridlines
rect(xleft=par("usr")[1],xright=par("usr")[2],
     ybottom=c(0,-1,-3:-4,-6:-7,-9:-10,-12:-16,-18:-21)-0.25,
     ytop=c(0,-1,-3:-4,-6:-7,-9:-10,-12:-16,-18:-21)+0.25,
     border=NA,col=cols.CTRL[4])
abline(v=log2(c(1,1.25,1.5,1.75,2,2.5,3)),col=cols.CTRL[3])
abline(v=log2(c(1,1.5,2)),col=cols.CTRL[2])
abline(v=0)

#####Add X-axis
axis(1,at=log2(c(1,1.25,1.5,1.75,2,2.5,3)),labels=NA)
axis(1,at=log2(c(1,1.25,1.5,1.75,2,2.5,3)),tick=F,line=-0.4,
     labels=c(1,1.25,1.5,1.75,2,2.5,3),cex.axis=0.75)

#####Add Y-Axis
axis(2,at=c(0,-1,-3:-4,-6:-7,-9:-10,-12:-16,-18:-21),labels=NA,tck=-0.015)

#####Plot points & 95% CIs
#Expected
segments(x0=EXP.mean-EXP.CI,x1=EXP.mean+EXP.CI,y0=0,y1=0,lwd=0.75)
points(x=EXP.mean,y=0,pch=21,bg=cols.CTRL,lwd=0.5,cex=1.3)
#Average across all gene sets
s=1
segments(x0=c(ALL.mean-ALL.CI,
              GERM.mean-GERM.CI,
              CNCR.mean-CNCR.CI),
         x1=c(ALL.mean+ALL.CI,
              GERM.mean+GERM.CI,
              CNCR.mean+CNCR.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.mean,GERM.mean,CNCR.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#OMIM association (average)
s=3
segments(x0=c(ALL.OMIM.avg.mean-ALL.OMIM.avg.CI,
              GERM.OMIM.avg.mean-GERM.OMIM.avg.CI,
              CNCR.OMIM.avg.mean-CNCR.OMIM.avg.CI),
         x1=c(ALL.OMIM.avg.mean+ALL.OMIM.avg.CI,
              GERM.OMIM.avg.mean+GERM.OMIM.avg.CI,
              CNCR.OMIM.avg.mean+CNCR.OMIM.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.OMIM.avg.mean,GERM.OMIM.avg.mean,CNCR.OMIM.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#OMIM (matched)
s=4
segments(x0=matched[[1]][c(1,3,5)]-matched[[1]][c(2,4,6)],
         x1=matched[[1]][c(1,3,5)]+matched[[1]][c(2,4,6)],
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=matched[[1]][c(1,3,5)],
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Specific expressors (average)
s=6
segments(x0=c(ALL.SpecExpr.avg.mean-ALL.SpecExpr.avg.CI,
              GERM.SpecExpr.avg.mean-GERM.SpecExpr.avg.CI,
              CNCR.SpecExpr.avg.mean-CNCR.SpecExpr.avg.CI),
         x1=c(ALL.SpecExpr.avg.mean+ALL.SpecExpr.avg.CI,
              GERM.SpecExpr.avg.mean+GERM.SpecExpr.avg.CI,
              CNCR.SpecExpr.avg.mean+CNCR.SpecExpr.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.SpecExpr.avg.mean,GERM.SpecExpr.avg.mean,CNCR.SpecExpr.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Specific expressors (matched)
s=7
segments(x0=matched[[4]][c(1,3,5)]-matched[[4]][c(2,4,6)],
         x1=matched[[4]][c(1,3,5)]+matched[[4]][c(2,4,6)],
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=matched[[4]][c(1,3,5)],
       y=-s+c(0.1,0,-0.1),
       pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#High expressors (average)
s=9
segments(x0=c(ALL.HighExpr.avg.mean-ALL.HighExpr.avg.CI,
              GERM.HighExpr.avg.mean-GERM.HighExpr.avg.CI,
              CNCR.HighExpr.avg.mean-CNCR.HighExpr.avg.CI),
         x1=c(ALL.HighExpr.avg.mean+ALL.HighExpr.avg.CI,
              GERM.HighExpr.avg.mean+GERM.HighExpr.avg.CI,
              CNCR.HighExpr.avg.mean+CNCR.HighExpr.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.HighExpr.avg.mean,GERM.HighExpr.avg.mean,CNCR.HighExpr.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#High expressors (matched)
s=10
segments(x0=matched[[2]][c(1,3,5)]-matched[[2]][c(2,4,6)],
         x1=matched[[2]][c(1,3,5)]+matched[[2]][c(2,4,6)],
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=matched[[2]][c(1,3,5)],
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Haplosufficient genes
s=12
segments(x0=c(ALL.Haplosuff.avg.mean-ALL.Haplosuff.avg.CI,
              GERM.Haplosuff.avg.mean-GERM.Haplosuff.avg.CI,
              CNCR.Haplosuff.avg.mean-CNCR.Haplosuff.avg.CI),
         x1=c(ALL.Haplosuff.avg.mean+ALL.Haplosuff.avg.CI,
              GERM.Haplosuff.avg.mean+GERM.Haplosuff.avg.CI,
              CNCR.Haplosuff.avg.mean+CNCR.Haplosuff.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.Haplosuff.avg.mean,GERM.Haplosuff.avg.mean,CNCR.Haplosuff.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Not haplosufficient genes
s=13
segments(x0=c(ALL.NotHaplosuff.avg.mean-ALL.NotHaplosuff.avg.CI,
              GERM.NotHaplosuff.avg.mean-GERM.NotHaplosuff.avg.CI,
              CNCR.NotHaplosuff.avg.mean-CNCR.NotHaplosuff.avg.CI),
         x1=c(ALL.NotHaplosuff.avg.mean+ALL.NotHaplosuff.avg.CI,
              GERM.NotHaplosuff.avg.mean+GERM.NotHaplosuff.avg.CI,
              CNCR.NotHaplosuff.avg.mean+CNCR.NotHaplosuff.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.NotHaplosuff.avg.mean,GERM.NotHaplosuff.avg.mean,CNCR.NotHaplosuff.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Constrained genes
s=14
segments(x0=c(ALL.Const.avg.mean-ALL.Const.avg.CI,
              GERM.Const.avg.mean-GERM.Const.avg.CI,
              CNCR.Const.avg.mean-CNCR.Const.avg.CI),
         x1=c(ALL.Const.avg.mean+ALL.Const.avg.CI,
              GERM.Const.avg.mean+GERM.Const.avg.CI,
              CNCR.Const.avg.mean+CNCR.Const.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.Const.avg.mean,GERM.Const.avg.mean,CNCR.Const.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Highly constrained genes
s=15
segments(x0=c(ALL.HiConst.avg.mean-ALL.HiConst.avg.CI,
              GERM.HiConst.avg.mean-GERM.HiConst.avg.CI,
              CNCR.HiConst.avg.mean-CNCR.HiConst.avg.CI),
         x1=c(ALL.HiConst.avg.mean+ALL.HiConst.avg.CI,
              GERM.HiConst.avg.mean+GERM.HiConst.avg.CI,
              CNCR.HiConst.avg.mean+CNCR.HiConst.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.HiConst.avg.mean,GERM.HiConst.avg.mean,CNCR.HiConst.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Extremely constrained genes
s=16
segments(x0=c(ALL.ExtConst.avg.mean-ALL.ExtConst.avg.CI,
              GERM.ExtConst.avg.mean-GERM.ExtConst.avg.CI,
              CNCR.ExtConst.avg.mean-CNCR.ExtConst.avg.CI),
         x1=c(ALL.ExtConst.avg.mean+ALL.ExtConst.avg.CI,
              GERM.ExtConst.avg.mean+GERM.ExtConst.avg.CI,
              CNCR.ExtConst.avg.mean+CNCR.ExtConst.avg.CI),
         y0=-s+c(0.1,0,-0.1),y1=-s+c(0.1,0,-0.1),lwd=0.75)
points(x=c(ALL.ExtConst.avg.mean,GERM.ExtConst.avg.mean,CNCR.ExtConst.avg.mean),
       y=-s+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
#Constratint-expression combos
s=18
sapply(1:4,function(i){
  segments(x0=c(combs.ALL.m[i]-combs.ALL.CI[i],
                combs.GERM.m[i]-combs.GERM.CI[i],
                combs.CNCR.m[i]-combs.CNCR.CI[i]),
           x1=c(combs.ALL.m[i]+combs.ALL.CI[i],
                combs.GERM.m[i]+combs.GERM.CI[i],
                combs.CNCR.m[i]+combs.CNCR.CI[i]),
           y0=-(s+i-1)+c(0.1,0,-0.1),y1=-(s+i-1)+c(0.1,0,-0.1),lwd=0.75)
  points(x=c(combs.ALL.m[i],combs.GERM.m[i],combs.CNCR.m[i]),
         y=-(s+i-1)+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]),lwd=0.5,cex=1.3)
})

#####Close device
dev.off()

# #Very high expressors (average)
# segments(x0=c(ALL.VeryExpr.avg.mean-ALL.VeryExpr.avg.CI,
#               GERM.VeryExpr.avg.mean-GERM.VeryExpr.avg.CI,
#               CNCR.VeryExpr.avg.mean-CNCR.VeryExpr.avg.CI),
#          x1=c(ALL.VeryExpr.avg.mean+ALL.VeryExpr.avg.CI,
#               GERM.VeryExpr.avg.mean+GERM.VeryExpr.avg.CI,
#               CNCR.VeryExpr.avg.mean+CNCR.VeryExpr.avg.CI),
#          y0=-4+c(0.1,0,-0.1),y1=-4+c(0.1,0,-0.1))
# points(x=c(ALL.VeryExpr.avg.mean,GERM.VeryExpr.avg.mean,CNCR.VeryExpr.avg.mean),
#        y=-4+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]))
# #Very high expressors (matched)
# segments(x0=matched[[3]][c(1,3,5)]-matched[[3]][c(2,4,6)],
#          x1=matched[[3]][c(1,3,5)]+matched[[3]][c(2,4,6)],
#          y0=-8+c(0.1,0,-0.1),y1=-8+c(0.1,0,-0.1))
# points(x=matched[[3]][c(1,3,5)],
#        y=-8+c(0.1,0,-0.1),pch=21,bg=c("white",cols.GERM[1],cols.CNCR[1]))
#Specific expressors (average)




