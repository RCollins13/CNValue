#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate correlation scatterplots for Fig3b-c

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

#####Load required libraries
require(MASS)





#################
#####CNCR vs GERM
#################
#####Read & clean data
GERM <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
GERM.means <- apply(GERM[,2:23],1,mean,na.rm=T)
GERM.means <- log2(GERM.means)
CNCR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
CNCR.means <- apply(CNCR[,-c(1:23)],1,mean,na.rm=T)
CNCR.means <- log2(CNCR.means)

#####Prepare plotting area
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/GERM_vs_CNCR.scatter.png",sep=""),
    width=4,height=4,units="in",res=1000)
par(mar=c(2,2,0.5,0.5))
plot(x=log2(c(1/5,9)),y=log2(c(1/5,9)),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")

#####Plot gridlines
abline(h=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),col=cols.CTRL[3],
       v=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)))
abline(h=0,v=0,lwd=1.5,col=cols.CTRL[1])

#####Plot x-axis
axis(1,at=log2(c(1/c(8:1),1:8)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(1,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),labels=NA)
axis(1,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),tick=F,
     labels=c("1/8","1/6","1/4","1/2",1,2,4,6,8),line=-0.4)

#####Plot y-axis
axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(2,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),labels=NA)
axis(2,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),tick=F,las=2,
     labels=c("1/8","1/6","1/4","1/2",1,2,4,6,8),line=-0.4)

#####Plot points
points(GERM.means,CNCR.means,
       pch=19,col=adjustcolor("black",alpha=0.5))

#####Plot density contours
means.df <- data.frame(GERM.means,CNCR.means)
means.df <- means.df[complete.cases(means.df),]
contours <- kde2d(means.df[,1],means.df[,2],n=200)
contour(contours,drawlabels=FALSE,nlevels=30,add=TRUE,
        col=adjustcolor(rev(rainbow(50)[1:40]),alpha=0.5))

#####Plot trendline
abline(lm(CNCR.means ~ GERM.means),lty=2)

#####Compute stats (add manually in illustrator)
cor(GERM.means,
    CNCR.means,
    use="complete.obs",
    method="spearman")
options(scipen=-1000)
cor.test(GERM.means,
         CNCR.means,
         method="spearman")$p.value
options(scipen=1000)

#####Close device
dev.off()







#################
#####DEL vs DUP
#################
#####Read & clean data
DEL <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "DEL_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DEL.means <- apply(DEL[,-1],1,mean,na.rm=T)
DEL.means <- log2(DEL.means)
DUP <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "DUP_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DUP.means <- apply(DUP[,-1],1,mean,na.rm=T)
mean(DUP.means,na.rm=T)
DUP.means <- log2(DUP.means)

#####Prepare plotting area
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/DEL_vs_DUP.scatter.png",sep=""),
    width=4,height=4,units="in",res=1000)
par(mar=c(2,2,0.5,0.5))
plot(x=log2(c(1/5,9)),y=log2(c(1/5,9)),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")

#####Plot gridlines
abline(h=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),col=cols.CTRL[3],
       v=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)))
abline(h=0,v=0,lwd=1.5,col=cols.CTRL[1])

#####Plot x-axis
axis(1,at=log2(c(1/c(8:1),1:8)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(1,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),labels=NA)
axis(1,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),tick=F,
     labels=c("1/8","1/6","1/4","1/2",1,2,4,6,8),line=-0.4)

#####Plot y-axis
axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
     labels=NA,tck=-0.01)
axis(2,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),labels=NA)
axis(2,at=log2(c(1/8,1/6,1/4,1/2,1,2,4,6,8)),tick=F,las=2,
     labels=c("1/8","1/6","1/4","1/2",1,2,4,6,8),line=-0.4)

#####Plot points
points(DEL.means,DUP.means,
       pch=19,col=adjustcolor("black",alpha=0.5))

#####Plot density contours
means.df <- data.frame(DEL.means,DUP.means)
means.df <- means.df[complete.cases(means.df),]
contours <- kde2d(means.df[,1],means.df[,2],n=200)
contour(contours,drawlabels=FALSE,nlevels=25,add=TRUE,
        col=adjustcolor(rev(rainbow(45)[1:32]),alpha=0.5))

#####Plot trendline
abline(lm(DUP.means ~ DEL.means),lty=2)

#####Compute stats (add manually in illustrator)
cor(DEL.means,
    DUP.means,
    use="complete.obs",
    method="spearman")
options(scipen=-1000)
cor.test(DEL.means,
         DUP.means,
         method="spearman")$p.value
options(scipen=1000)

#####Close device
dev.off()





