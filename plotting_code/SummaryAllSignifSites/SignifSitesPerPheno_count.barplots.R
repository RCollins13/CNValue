#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to plot counts of significant sites by class per phenotype by DEL/DUP

#####Set parameters & load libraries
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")

#####Read & clean data
x <- read.table(paste(WRKDIR,"plot_data/SummaryAllSignifSites/signifSites.countPerPheno.txt",sep=""),
                header=F,sep="\t")
rownames(x) <- x[,1]
x <- t(apply(x[,-1],1,as.integer))
colnames(x) <- c("DEL.seg","DEL.gene","DEL.reg",
                 "DUP.seg","DUP.gene","DUP.reg")
x <- x[nrow(x):1,]

#####PLOT
#Open device
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/SummaryAllSignifSites/",
          "SignifSitesPerPheno_count.barplots.pdf",sep=""),
    height=2.5,width=5)
#Set plot parameters
xmax <- 330
alphas <- c(1,0.6,0.2)
DEL.cols <- sapply(alphas,function(a){
  return(adjustcolor("red",alpha=a))
})
DUP.cols <- sapply(alphas,function(a){
  return(adjustcolor("blue",alpha=a))
})
par(mar=c(0.5,4.5,2.5,1.5),bty="n")
plot(x=c(0,xmax),y=c(0,nrow(x)),type="n",
     xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

#Add gridlines
abline(v=seq(0,xmax,15),lwd=0.7,col=cols.CTRL[3])
abline(v=seq(0,xmax,30),col=cols.CTRL[1])
abline(v=0,lwd=3)

#Print axes
axis(2,at=(1:(nrow(x)-1))-0.5,tick=F,labels=rownames(x)[1:(nrow(x)-1)],las=2,line=0.3)
axis(2,at=nrow(x)-0.5,tick=F,labels=rownames(x)[nrow(x)],las=2,line=-0.9)
axis(3,at=seq(0,xmax,15),labels=NA,tck=-0.01,col=cols.CTRL[1])
axis(3,at=seq(0,xmax,30),labels=NA,tck=-0.025)
axis(3,at=seq(30,xmax,60),tick=F,line=-0.8)
axis(3,at=seq(0,xmax,60),tick=F,line=-0.8)
mtext(3,line=1.2,text="Significant Loci (Count)")

#Plot DEL and DUP rectangles
sapply(1:nrow(x),function(i){
  #Get counts
  DEL.counts <- as.numeric(c(0,cumsum(x[i,1:3])))
  DUP.counts <- as.numeric(c(0,cumsum(x[i,4:6])))

  #Plot rectangles
  rect(xleft=0,xright=max(DEL.counts),
       ybottom=i-0.5,ytop=i-0.2,col="white")
  rect(xleft=DEL.counts[1:3],xright=DEL.counts[2:4],
       ybottom=i-0.5,ytop=i-0.2,col=DEL.cols)
  rect(xleft=0,xright=max(DUP.counts),
       ybottom=i-0.8,ytop=i-0.5,col="white")
  rect(xleft=DUP.counts[1:3],xright=DUP.counts[2:4],
       ybottom=i-0.8,ytop=i-0.5,col=DUP.cols)
})

#Close device
dev.off()




