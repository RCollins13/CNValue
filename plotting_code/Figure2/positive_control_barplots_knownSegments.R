#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate barplots of large segment positive control sets

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
# cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

######################################################
#####Read data, clean, and calculate binomial p-values
######################################################
dat <- read.table(paste(WRKDIR,"plot_data/figure2/known_locus_overlap.txt",sep=""),header=T)
plotDat <- t(sapply(seq(3,ncol(dat),2),function(i){
  matrix(unlist(sapply(1:2,function(r){
    vals <- as.numeric(dat[r,c(2,i,i+1)])
    vals <- c(vals,vals[2]/vals[1],vals[3]/vals[1],vals[2]/vals[3])
    binom.res <- binom.test(x=vals[2],n=vals[1],p=vals[5],alternative="greater")
    vals <- c(vals,binom.res$conf.int[1:2],binom.res$p.value)
    return(vals)
  })),nrow=1)
}))
rownames(plotDat) <- unique(gsub(".expected","",gsub(".observed","",colnames(dat[,-c(1:2)]))))
colnames(plotDat) <- as.vector(sapply(c("DEL","DUP"),function(CNV){
  return(paste(CNV,c("total","obs","exp","obsFrac","expFrac","fold","lower","upper","p"),sep="."))
}))
#Subset to values for plotting
plotDat <- as.data.frame(plotDat[rev(c(1,3,13,4,6,7,9,2)),])
# plotDat <- plotDat[rev(c(1,1,3,13,4,6,7,9,2)),]
# plotDat[nrow(plotDat)-1,] <- NA
# rownames(plotDat)[nrow(plotDat)-1] <- NA
#Set plot labels
# plotLab <- c("Marshall (2016)","Schaefer (2013)",
#              "Pinto (2014)","DDD","Coe (2014)",
#              "Wapner (2012)","ClinGen","","All studies")
plotLab <- c("Marshall (2016)","Schaefer (2013)",
             "Pinto (2014)","DDD","Coe (2014)",
             "Wapner (2012)","ClinGen","All studies")

#################
######Plot values
#################
buffer=0.2
#Prep plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/largeSegments_positiveControls.barplot.pdf",sep=""),
    height=4.5,width=5)
par(mar=c(0.5,7,3,1),bty="n")
plot(x=c(0,2+buffer),y=c(0,nrow(plotDat)),type="n",
     xaxs="i",xlab="",xaxt="n",yaxs="i",ylab="",yaxt="n")

#Iterate over data and plotprint stats
sapply(1:nrow(plotDat),function(i){
  if(!any(is.na(plotDat[i,]))){
    rect(xleft=c(0,0,1+buffer,1+buffer),
         xright=c(0,0,1+buffer,1+buffer)+1,
         ybottom=rep(c(i-0.9,i-0.5),2),
         ytop=rep(c(i-0.5,i-0.1),2),
         col="#F1F1F2",border="#A5A6A7")
    segments(x0=c(seq(0.05,0.95,0.05),seq(1.05,1.95,0.05)+buffer),
             x1=c(seq(0.05,0.95,0.05),seq(1.05,1.95,0.05)+buffer),
             y0=i-0.9,y1=i-0.1,col="#E3E4E5")
    rect(xleft=c(0,0,1+buffer,1+buffer),
         xright=c(0,0,1+buffer,1+buffer)+1,
         ybottom=rep(c(i-0.9,i-0.5),2),
         ytop=rep(c(i-0.5,i-0.1),2),
         col=NA,border="#A5A6A7")
    rect(xleft=c(0,0,1+buffer,1+buffer),
         xright=c(plotDat[i,c(4,5)],plotDat[i,c(13,14)]+1+buffer),
         ybottom=c(i-0.5,i-0.9,i-0.5,i-0.9),
         ytop=c(i-0.1,i-0.5,i-0.1,i-0.5),
         col=c("red","#373738","blue","#373738"))
    segments(x0=seq(0.05,plotDat[i,4],0.05),x1=seq(0.05,plotDat[i,4],0.05),
             y0=i-0.5,y1=i-0.1,col="#cc0000")
    segments(x0=seq(0.05,plotDat[i,13],0.05)+1+buffer,
             x1=seq(0.05,plotDat[i,13],0.05)+1+buffer,
             y0=i-0.5,y1=i-0.1,col="#000099")
    if(plotDat[i,5]>0.05){
      segments(x0=seq(0.05,plotDat[i,5],0.05),x1=seq(0.05,plotDat[i,5],0.05),
               y0=i-0.9,y1=i-0.5)
    }
    if(plotDat[i,14]>0.05){
      segments(x0=seq(0.05,plotDat[i,14],0.05)+1+buffer,
               x1=seq(0.05,plotDat[i,14],0.05)+1+buffer,
               y0=i-0.9,y1=i-0.5)
    }
    rect(xleft=c(0,0,1+buffer,1+buffer),
         xright=c(plotDat[i,c(4,5)],plotDat[i,c(13,14)]+1+buffer),
         ybottom=c(i-0.5,i-0.9,i-0.5,i-0.9),
         ytop=c(i-0.1,i-0.5,i-0.1,i-0.5),col=NA)
    # text(x=c(1,1,2+buffer,2+buffer),
    #      y=c(i-0.3,i-0.7,i-0.3,i-0.7),
    #      labels=c(paste(plotDat[i,2],"/",plotDat[i,1],sep=""),
    #               paste(plotDat[i,3],"/",plotDat[i,1],sep=""),
    #               paste(plotDat[i,11],"/",plotDat[i,10],sep=""),
    #               paste(plotDat[i,12],"/",plotDat[i,10],sep="")),
    #      col=c("red","#6E6F70","blue","#6E6F70"),pos=4,cex=0.8)
  }
})

#Add study labels
axis(2,at=c(1:nrow(plotDat))-0.5,tick=F,line=-0.8,labels=plotLab,las=2,cex.axis=1.2)

#Add x-axis
axis(3,at=seq(0,1,0.05),labels=NA,tck=-0.01,col="#A5A6A7")
axis(3,at=seq(1+buffer,2+buffer,0.05),labels=NA,tck=-0.01,col="#A5A6A7")
axis(3,at=seq(0,1,0.2),labels=NA,tck=-0.02)
axis(3,at=seq(1,2,0.2)+buffer,labels=NA,tck=-0.02)
axis(3,at=c(0,0.4,0.8,1+buffer,1.4+buffer,1.8+buffer),
     labels=format(c(0,0.4,0.8,0,0.4,0.8),nsmall=1),
     tick=F,line=-0.7,cex.axis=0.8)
axis(3,at=c(0.2,0.6,1,1.2+buffer,1.6+buffer,2+buffer),
     labels=format(c(0.2,0.6,1,0.2,0.6,1),nsmall=1),
     tick=F,line=-0.7,cex.axis=0.8)

#Close device
dev.off()



