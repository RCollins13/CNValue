#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate reverse CDFs for CNV size

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
logvect.all <- log10(as.vector(sapply(4:7,function(x){return((1:9)*10^x)})))
logvect.mids <- log10(as.vector(sapply(4:6,function(x){return(c(5,10)*10^x)})))
logvect.mains <- log10(as.vector(sapply(4:6,function(x){return(10^x)})))
pctvect.all <- log10(as.vector(sapply(-1:2,function(x){return((1:9)*10^x)})))
pctvect.mids <- log10(as.vector(sapply(-1:1,function(x){return(c(5,10)*10^x)})))
pctvect.mains <- log10(as.vector(sapply(-1:2,function(x){return(10^x)})))
cols.groups <- c("#A5A6A7","#00BFF4","#EC008D","#FFCB00")
phenos <- c("NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
cols.phenos <- c(rep(cols.groups[2],10),
                 rep(cols.groups[3],10),
                 rep(cols.groups[4],13))

#####Master function
plotsizes <- function(filt,CNV,cols,xaxis=T,yaxis=T,subgroups=T,
                      mar=c(3.5,3.5,0.5,0.5)){
  #Read data & convert to log-scaled CDF
  sizes <- lapply(list("CTRL","NEURO","SOMA","CNCR"),function(group){
    dat <- read.table(paste(WRKDIR,"plot_data/figure1/",CNV,"_size.",group,
                            ".noMaxSize.E2.",filt,".txt",sep=""))[,1]
    dat <- log10(dat)
    cdf <- sapply(logvect.all,function(min){
      length(which(dat>=min))/length(dat)
    })
    cdf <- log10(100*cdf)
    cdf[which(!is.finite(cdf))] <- -100
    return(cdf)
  })

  #Read data & compute medians/IQRs
  distribs <- lapply(list("CTRL","NEURO","SOMA","CNCR"),function(group){
    dat <- read.table(paste(WRKDIR,"plot_data/figure1/",CNV,"_size.",group,
                            ".noMaxSize.E2.",filt,".txt",sep=""))[,1]
    dat <- log10(dat)

    return(as.vector(summary(dat)[c(2:3,5)]))
  })
  medians <- unlist(lapply(distribs,function(v){return(v[2])}))
  IQRs <- t(as.data.frame(lapply(distribs,function(v){return(v[c(1,3)])})))

  #Read subgroup data if optioned
  if(subgroups==T){
    #Load data & compute CDFs
    sizes.sub <- lapply(phenos,function(group){
      dat <- read.table(paste(WRKDIR,"plot_data/figure1/",CNV,"_size.",group,
                              ".noMaxSize.E2.",filt,".txt",sep=""))[,1]
      dat <- log10(dat)
      cdf <- sapply(logvect.all,function(min){
        length(which(dat>=min))/length(dat)
      })
      cdf <- log10(100*cdf)
      cdf[which(!is.finite(cdf))] <- -100
      return(cdf)
    })
  }

  #Set parameters
  par(mar=mar,bty="n")

  #Prepare plot
  plot(x=log10(c(50000,5000000)),y=log10(c(0.03,100)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Gridlines - vertical & horizontal interleaved
  abline(v=logvect.all,col="gray95")
  abline(h=pctvect.all,col="gray95")
  abline(v=logvect.mids,col="gray85")
  abline(h=pctvect.mids,col="gray85")
  abline(v=logvect.mains,col="gray75")
  abline(h=pctvect.mains,col="gray75")

  #X-axis
  axis(1,logvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
  axis(1,logvect.mids,labels=NA,tck=-0.02,col="gray20",lwd=1)
  axis(1,logvect.mains,labels=NA,tck=-0.025,col="black",lwd=1.5)
  if(xaxis==T){
    axis(1,at=logvect.mids,line=-0.6,tick=F,las=2,cex.axis=0.7,
         labels=c("50kb","100kb","500kb","1Mb","5Mb","10Mb"))
  }

  #Y-axis
  axis(2,pctvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
  axis(2,pctvect.mids,labels=NA,tck=-0.02,col="gray20",lwd=1)
  axis(2,pctvect.mains,labels=NA,tck=-0.025,col="black",lwd=1.5)
  if(yaxis==T){
    axis(2,at=c(log10(0.1),pctvect.mids),line=-0.5,tick=F,las=2,
         labels=paste(c(0.1,0.5,1,5,10,50,100),"%",sep=""))
  }

  #Plot CDFs for subgroups, if optioned
  if(subgroups==T){
    sapply(1:length(sizes.sub),function(i){
      points(logvect.all,sizes.sub[[i]],type="l",col=cols.phenos[i])
    })
  }

  #Plot CDFs
  sapply(1:length(sizes),function(i){
    points(logvect.all,sizes[[i]],type="l",lwd=4)
    points(logvect.all,sizes[[i]],type="l",col=cols[i],lwd=3)
  })

  #Medians & IQRs
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1),
       col="white",border="white")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1),ytop=par("usr")[4])
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1)-0.1,col="white")
  segments(x0=logvect.all,x1=logvect.all,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray95")
  segments(x0=logvect.mids,x1=logvect.mids,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray85")
  segments(x0=logvect.mains,x1=logvect.mains,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray75")
  rect(xleft=rev(IQRs[,1]),xright=rev(IQRs[,2]),
       ybottom=seq(par("usr")[3],log10(0.1)-0.1,
                   by=(log10(0.1)-par("usr")[3])/5)[1:4],
       ytop=seq(par("usr")[3],log10(0.1)-0.1,
                by=(log10(0.1)-par("usr")[3])/5)[2:5],
       border=NA,col=rev(adjustcolor(cols,alpha=0.75)))
  # points(x=medians,y=rep(mean(c(par("usr")[3],log10(0.1)-0.1)),4),
  #        pch=18,col=cols,cex=2)
  median.pos <- sapply(1:4,function(i){
    pts <- rev(seq(par("usr")[3],log10(0.1)-0.1,
                   by=(log10(0.1)-par("usr")[3])/5))
    mean(pts[c(i,i+1)])
  })
  points(x=medians,y=median.pos,pch=23,cex=0.5,col="black",bg=cols,lwd=0.6)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1)-0.1,ytop=log10(0.1),
       col="white",border="white")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1)-0.1,col=NA)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1),ytop=par("usr")[4],col=NA)
}

#Run plots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/CNVsize_by_filter_by_group.reverse_CDF.pdf",sep=""),
    width=7,height=1.8)
par(mfrow=c(1,4))
plotsizes(filt="coding",CNV="DEL",cols=cols.groups,subgroups=T)
plotsizes(filt="coding",CNV="DUP",cols=cols.groups,subgroups=T)
plotsizes(filt="noncoding",CNV="DEL",cols=cols.groups,subgroups=T)
plotsizes(filt="noncoding",CNV="DUP",cols=cols.groups,subgroups=T)
dev.off()



