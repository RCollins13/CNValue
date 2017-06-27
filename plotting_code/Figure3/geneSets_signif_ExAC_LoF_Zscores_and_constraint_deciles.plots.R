#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate ExAC LoF z-scores for genes in significant gene sets for Fig3

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
cols.CNV <- c("gray20","red","blue")

#####Load libraries
require(beeswarm)
require(plotrix)
require(vioplot)

#####Read data
dat <- lapply(list("allGenes","noPhenos","anyPhenos"),function(category){
  dat <- read.table(paste(WRKDIR,"plot_data/figure3/ExAC_LoF_z.",category,".txt",sep=""),
                    header=F)[,1]
  return(dat)
})

#####Prepare plotting area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/geneSet_ExAC_z_scores.boxplots.pdf",sep=""),
    width=1.5,height=3.5)
par(bty="n",mar=c(1,2.5,0.5,0.5))
plot(x=c(0.5,3.5),y=c(-3,7),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")

#####Add boxplots
boxplot(dat,add=T,notch=T,outline=F,lty=1,staplewex=0,
        xaxt="n",yaxt="n",col=cols.CTRL[2])

#####Add axes
axis(1,at=1:3,labels=NA)
axis(2,at=-6:8,labels=NA,tck=-0.013,col=cols.CTRL[1])
axis(2,at=seq(-6,8,2),labels=NA)
axis(2,at=seq(-6,8,2),tick=F,las=2,line=-0.3,cex.axis=1.25)
# mtext(2,text="LoF Constraint Z-Score",line=1.3)

#Close device
dev.off()

#####Run t-tests between groups
options(scipen=-100)
t.test(dat[[1]],dat[[3]],alternative="less")$p.value
options(scipen=1000)




#################################################
#####Dotplots of rCNV ORs at constraint quintiles
#################################################

#####Read data
ORs <- lapply(c("CNV","DEL","DUP"),function(CNV){
  OR <- read.table(paste(WRKDIR,"plot_data/figure3/constraint_quintiles/",
                         CNV,"_E2_all_exonic.effectSizes.txt",sep=""),header=F)
  names(OR) <- c("anno",phenos)
  return(OR)
})

#####Split by germline & cancer
GERM <- lapply(ORs,function(df){
  return(df[,c(1:23)])
})
CNCR <- lapply(ORs,function(df){
  return(df[,-c(1:23)])
})

#####Master function to create dotplot
plotOR <- function(ORs){
  #####Get mean & 95% CI for each decile
  stats <- lapply(ORs,function(df){
    dat <- apply(df[,-1],1,function(row){
      m <- mean(log2(row),na.rm=T)
      ci <- 1.96*std.error(log2(row),na.rm=T)
      return(c(m,m-ci,m+ci))
    })
    return(as.data.frame(t(dat)))
  })

  #####Prepare plotting area
  par(mar=c(2.5,2.5,0.5,1))
  plot(x=c(1,5),y=log2(c(2/3,2.5)),type="n",
       yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #####Draw gridlines
  abline(h=log2(c(1/3,1/2,3/4,1,1.5,2,2.5,3)),col=cols.CTRL[2])
  abline(h=0,lwd=2)

  #####Iterate over stats and plot CNV/DEL/DUP
  pad <- c(-0.05,0,0.05)
  sapply(1:3,function(i){
    #Trendlines
    segments(x0=(1:4)+pad[i],x1=(2:5)+pad[i],
             y0=stats[[i]][1:4,1],y1=stats[[i]][2:5,1],
             col=cols.CNV[i])
    #95% Confidence Intervals
    segments(y0=stats[[i]][,2],y1=stats[[i]][,3],
             x0=(1:5)+pad[i],x1=(1:5)+pad[i],
             col=cols.CNV[i])
    #Point estimates
    points(y=stats[[i]][,1],x=(1:5)+pad[i],
           pch=21,bg=cols.CNV[i])
  })

  #####Add Y axis
  axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
       labels=NA,tck=-0.01)
  axis(2,at=log2(c(1/6,1/4,1/2,3/4,1,1.5,2,2.5,3,4,6)),labels=NA)
  axis(2,at=log2(c(1/6,1/4,1/2,3/4,1,1.5,2,2.5,3,4,6)),tick=F,line=-0.3,
       labels=c("1/6","1/4","1/2","3/4",1,1.5,2,2.5,3,4,6),las=2)

  #####Add X-Axis
  # axis(1,at=1:5,labels=NA)
  axis(1,at=1:5,tick=F,line=-1,labels=c("0-20","20-40","40-60","60-80","80-100"),cex.axis=0.75)
}


#####Plot ORs
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/ExAC_constraint_quintiles.GERM.pdf",sep=""),
    width=3,height=2.3)
plotOR(GERM)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/ExAC_constraint_quintiles.CNCR.pdf",sep=""),
    width=3,height=2.3)
plotOR(CNCR)
dev.off()



###############################################
#####Dotplots of rCNV ORs at constraint deciles
###############################################

#####Read data
ORs <- lapply(c("CNV","DEL","DUP"),function(CNV){
  OR <- read.table(paste(WRKDIR,"plot_data/figure3/constraint_deciles/",
                         CNV,"_E2_all_exonic.effectSizes.txt",sep=""),header=F)
  names(OR) <- c("anno",phenos)
  return(OR)
})

#####Split by germline & cancer
GERM <- lapply(ORs,function(df){
  return(df[,c(1:23)])
})
CNCR <- lapply(ORs,function(df){
  return(df[,-c(1:23)])
})

#####Master function to create dotplot
plotOR.deciles <- function(ORs){
  #####Get mean & 95% CI for each decile
  stats <- lapply(ORs,function(df){
    dat <- apply(df[,-1],1,function(row){
      m <- mean(log2(row),na.rm=T)
      ci <- 1.96*std.error(log2(row),na.rm=T)
      return(c(m,m-ci,m+ci))
    })
    return(as.data.frame(t(dat)))
  })

  #####Prepare plotting area
  par(mar=c(2.5,2.5,0.5,1))
  plot(x=c(-1,-10),y=log2(c(0.5,3)),type="n",
       yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #####Draw gridlines
  abline(h=log2(c(1/3,1/2,3/4,1,1.5,2,2.5,3)),col=cols.CTRL[2])
  abline(h=0)

  #####Iterate over stats and plot CNV/DEL/DUP
  pad <- c(-0.05,0,0.05)
  sapply(2:3,function(i){
    #Trendlines
    linfit.dat <- data.frame(-1:-10,stats[[i]][,1])
    names(linfit.dat) <- c("decile","estimate")
    abline(lm(estimate ~ decile, data=linfit.dat),
           col=cols.CNV[i],lwd=1.8)
    #Connecting lines
    points(y=stats[[i]][,1],x=(-1:-10)-pad[i],
           type="l",lty=5,col=cols.CNV[i],lwd=0.85)
  })
  sapply(2:3,function(i){
    #95% Confidence Intervals
    segments(y0=stats[[i]][,2],y1=stats[[i]][,3],
             x0=(-1:-10)-pad[i],x1=(-1:-10)-pad[i])
    #Point estimates
    points(y=stats[[i]][,1],x=(-1:-10)-pad[i],
           pch=21,bg=cols.CNV[i],lwd=0.5)
  })

  #####Add Y axis
  axis(2,at=log2(c(1/c(6:1),1:6)),col=cols.CTRL[1],
       labels=NA,tck=-0.01)
  axis(2,at=log2(c(1/6,1/4,1/2,3/4,1,1.5,2,2.5,3,4,6)),labels=NA)
  axis(2,at=log2(c(1/6,1/4,1/2,3/4,1,1.5,2,2.5,3,4,6)),tick=F,line=-0.3,
       labels=c("1/6","1/4","1/2","3/4",1,1.5,2,2.5,3,4,6),las=2)

  #####Add X-Axis
  # axis(1,at=1:5,labels=NA)
  # axis(1,at=-1:-10,tick=F,line=-0.9,cex.axis=0.75,las=2,
  #      labels=rev(c("1st","2nd","3rd",paste(4:10,"th",sep=""))))
}


#####Plot ORs
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/ExAC_constraint_deciles.GERM.pdf",sep=""),
    width=3,height=2.3)
plotOR.deciles(GERM)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/ExAC_constraint_deciles.CNCR.pdf",sep=""),
    width=3,height=2.3)
plotOR.deciles(CNCR)
dev.off()






