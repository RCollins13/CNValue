#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate violin plots of CNV size per VF filter

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
logvect.all <- log10(as.vector(sapply(4:7,function(x){return((1:9)*10^x)})))
logvect.mids <- log10(as.vector(sapply(4:6,function(x){return(c(5,10)*10^x)})))
logvect.mains <- log10(as.vector(sapply(4:7,function(x){return(10^x)})))
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(vioplot)

#####Master function
plotsizes <- function(group,cols,xaxis=T,yaxis=T,h=1/80,line=T,
                      mar=c(2.5,2.5,0.5,0.5)){
  #Read data & convert to log-scaled CDF
  sizes <- lapply(list("E2","E3","E4","N1"),function(freq){
    sizes <- lapply(list("DEL","DUP"),function(CNV){
      dat <- read.table(paste(WRKDIR,"plot_data/figure1/",CNV,"_size.",group,
                              ".",freq,".all.txt",sep=""))[,1]
      return(log10(dat))
    })
    return(sizes)
  })

  #Set parameters
  par(mar=mar,bty="n")

  #Prepare plot
  plot(x=c(0.5,4.5),y=log10(c(50000,5000000)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Gridlines - horizontal
  abline(h=logvect.all,col="gray95")
  abline(h=logvect.mids,col="gray85")
  abline(h=logvect.mains,col="gray75")

  #Y-axis
  axis(2,logvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
  axis(2,logvect.mids,labels=NA,tck=-0.02,col="gray20",lwd=1)
  axis(2,logvect.mains,labels=NA,tck=-0.025,col="black",lwd=1.5)
  if(yaxis==T){
    axis(2,at=logvect.mids,line=-0.6,tick=F,las=2,cex.axis=0.9,
         labels=c("50kb","100kb","500kb","1Mb","5Mb","10Mb"))
  }

  #X-axis
  if(xaxis==T){
    axis(1,at=1:4,line=-1,tick=F,labels=c("1%","0.1%","0.01%","n=1"),cex.axis=0.8)
    mtext(1,line=0.8,text="Maximum VF",cex=0.6)
  }

  #Plot violins
  stats <- sapply(1:length(sizes),function(i){
    vioplot(unlist(sizes[[i]]),h=h,names=NA,col=cols[i],at=i,add=T,pchMed=18)
  })

  #Plot borders around medians
  points(x=1:4,y=unlist(stats[3,]),pch=23,col="black",bg="white",cex=1.2)

  #Fit line
  if(line==T){
    stats <- as.data.frame(t(stats))
    stats <- data.frame(unlist(stats$median),1:4)
    stats$point <- 1:4
    names(stats) <- c("median","VF")
    abline(lm(median ~ VF, stats))
    summary(lm(median ~ VF, stats))
  }

}

#Run plots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/CNVsize_by_VF.violins.pdf",sep=""),
    width=8,height=1.6)
par(mfrow=c(1,4))
plotsizes(group="CTRL",cols=cols.CTRL)
plotsizes(group="NEURO",cols=cols.NEURO)
plotsizes(group="SOMA",cols=cols.SOMA)
plotsizes(group="CNCR",cols=cols.CNCR,xaxis=F,h=1/100,line=F)
dev.off()







#####Code to plot bars for del/dup per person per bin
plotbars <- function(group,n,xaxis=T,mar=c(0.5,2.5,1.5,0.5)){
  #Read data & convert to log-scaled CDF
  pp <- lapply(list("E2","E3","E4","N1"),function(freq){
    pp <- lapply(list("DEL","DUP"),function(CNV){
      dat <- read.table(paste(WRKDIR,"plot_data/figure1/",CNV,"_size.",group,
                              ".",freq,".all.txt",sep=""))[,1]
      return(length(dat)/n)
    })
    return(unlist(pp))
  })

  #Set parameters
  par(mar=mar,bty="n")

  #Prepare plot
  plot(x=c(0.5,4.5),y=c(0,max(unlist(pp))),
       type="n",xaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",las=2)

  #Plot bars
  sapply(1:length(pp),function(i){
    rect(xleft=c(i-0.225,i),xright=c(i,i+0.225),
         ybottom=0,ytop=pp[[i]],col=c("red","blue"),border=NA)
  })

  #Add x-axis
  if(xaxis==T){
    axis(3,at=1:4,line=-1,tick=F,labels=c("1%","0.1%","0.01%","n=1"),cex.axis=0.8)
  }

  #Add y-axis
  axis(2,at=axTicks(2)[c(1,3,5)],tck=-0.025,labels=F)
  axis(2,at=axTicks(2)[c(1,3,5)],tick=F,line=-0.7,las=2,cex.axis=0.7)
}

#Run plots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/DEL_DUP_split_by_VF.barplots.pdf",sep=""),
    width=8,height=0.7)
par(mfrow=c(1,4))
plotbars(group="CTRL",n=38628)
plotbars(group="NEURO",n=57760)
plotbars(group="SOMA",n=13361)
plotbars(group="CNCR",n=10844,xaxis=F)
dev.off()

