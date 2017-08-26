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
BOTH <- 11
GERM <- 15-BOTH
CNCR <- 53-BOTH

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

##################
#####Plotting code
##################
#Open device
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_GERM_VS_CNCR.barplot.pdf",sep=""),
    height=3,width=1.6)
#Prepare plot area
par(mar=c(0.5,3.1,0.5,0.1),bty="n")
plot(x=c(0,1),y=c(0,sum(BOTH,GERM,CNCR)),
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

#Draw gridlines
abline(h=seq(0,100,5),col=cols.CTRL[4])
abline(h=seq(0,100,10),col=cols.CTRL[3])

#Draw rectangles
rect(xleft=0.2,xright=0.8,
     ybottom=par("usr")[3],ytop=CNCR,
     col=cols.CNCR[1])
rect(xleft=0.2,xright=0.8,
     ybottom=CNCR,ytop=CNCR+BOTH,
     col=cols.GERM[1])
rect(xleft=0.2,xright=0.8,
     ybottom=CNCR,ytop=CNCR+BOTH,
     col=cols.CNCR[1],density=8,lwd=6,border=NA)
rect(xleft=0.2,xright=0.8,
     ybottom=CNCR,ytop=CNCR+BOTH,
     col=NA)
rect(xleft=0.2,xright=0.8,
     ybottom=CNCR+BOTH,ytop=CNCR+BOTH+GERM,
     col=cols.GERM[1])

#Clean up rectangles and draw gridlines
rect(xleft=c(par("usr")[1],0.8),
     xright=c(0.2,par("usr")[2]),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     border=NA,col="white")
segments(x0=c(rep(par("usr")[1],times=length(seq(0,100,5))),
              rep(0.8,times=length(seq(0,100,5))))
         x1=c(rep(0.2,times=length(seq(0,100,5))),
              rep(par("usr")[2],times=length(seq(0,100,5)))),
         y0=rep(seq(0,100,5),2),y1=rep(seq(0,100,5),2),
         col=cols.CTRL[4])

#Add axes
axis(1,at=0:1,labels=NA,tck=0)
axis(2,at=seq(0,100,5),col=cols.CTRL[1],tck=-0.02,labels=NA)
axis(2,at=seq(0,100,10),tck=-0.04,labels=NA)
axis(2,at=seq(0,100,10),tick=F,las=2,line=-0.6,cex.axis=1.1)

#Add cleanup box
rect(xleft=0.2,xright=0.8,
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col=NA)

#Close device
dev.off()
