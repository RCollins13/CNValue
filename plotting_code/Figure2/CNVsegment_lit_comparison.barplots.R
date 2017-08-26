#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to make barplots of overlap between significant rCNV segments and previous reports

####################################
#####Set parameters & load libraries
####################################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)

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
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_vs_knownSites.barplots.pdf",sep=""),
    height=3.6,width=2.5)
#Prepare plot area
par(mar=c(0.5,3.1,0.5,0.1),bty="n")
plot(x=c(0,2),y=c(0,1),
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

#Draw gridlines
abline(h=seq(0,1,0.1),col=cols.CTRL[4])
abline(h=seq(0,1,0.2),col=cols.CTRL[3])

#Add axes
axis(1,at=0:2,labels=NA)
axis(2,at=seq(0,1,0.1),col=cols.CTRL[1],tck=-0.01,labels=NA)
axis(2,at=seq(0,1,0.2),tck=-0.02,labels=NA)
axis(2,at=seq(0,1,0.2),tick=F,las=2,line=-0.6,cex.axis=1.1,
     labels=paste(seq(0,100,20),"%",sep=""))

#Draw rectangles
rect(xleft=c(0.2,1.2),xright=c(0.8,1.8),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col=cols.CTRL[1])
rect(xleft=c(0.2,1.2),xright=c(0.8,1.8),
     ybottom=par("usr")[3],ytop=c(12/15,38/53),
     col=c(cols.GERM[1],cols.CNCR[1]))

#Close device
dev.off()
