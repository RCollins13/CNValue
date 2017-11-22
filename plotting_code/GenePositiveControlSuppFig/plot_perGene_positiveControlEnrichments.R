#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate positive control barplots for per-gene burden analyses

#####Set parameters & load libraries
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)

#####Read data & compute statistics
x <- read.table(paste(WRKDIR,"plot_data/GenePositiveControlSuppFig/signifGenes_positiveControlData.txt",sep=""),header=F)
stats <- t(apply(x[,-1],1,function(vals){
  vals <- as.numeric(vals)
  fracs <- c(vals[2]/vals[1],
             vals[4]/vals[3],
             vals[6]/vals[5])
  folds <- c(fracs[2]/fracs[1],
             fracs[3]/fracs[1])
  if(!is.na(vals[2])){
    fold.CIs <- c(as.numeric(binom.test(vals[4],vals[3],p=vals[2]/vals[1])$conf.int/(vals[2]/vals[1])),
                  as.numeric(binom.test(vals[6],vals[5],p=vals[2]/vals[1])$conf.int/(vals[2]/vals[1])))
  }else{
    fold.CIs <- rep(NA,4)
  }
  p <- c(phyper(vals[4],vals[2],vals[1]-vals[2],vals[3],lower.tail=F),
         phyper(vals[6],vals[2],vals[1]-vals[2],vals[5],lower.tail=F))
  return(c(vals,fracs,folds,fold.CIs,p))
}))
x <- as.data.frame(cbind(x[,1],stats))
names(x) <- c("set","n.all","n.all_in_set","n.DEL","n.DEL_in_set",
              "n.DUP","n.DUP_in_set","frac.all","frac.DEL","frac.DUP",
              "fold.DEL","fold.DUP","lower.DEL","upper.DEL",
              "lower.DUP","upper.DUP","p.DEL","p.DUP")
x <- x[nrow(x):1,]

#####MASTER PLOT
#Set plot params
fold.max <- 10
set.labels <- c("LoF Constrained (ExAC; pLI>0.9)","LoF Tolerant (ExAC; pLI<0.1)",
                "","Haploinsufficient (ClinGen)","Disease Associated (ClinVar)",
                "Recurrent dnLoF in DDs (DDD)","DD Associated (Nguyen; extTADA)",
                "ASD Associated (Nguyen; extTADA)","ID Associated (Nguyen; extTADA)",
                "","Autosomal Dominant Disease (Berg + Blekhman)","Dominant LoF in DDs (DECIPHER)",
                "Dominant GoF in DDs (DECIPHER)","","Autosomal Recessive Disease (Berg + Blekhman)",
                "Recessive LoF in DDs (DECIPHER)","Recessive GoF in DDs (DECIPHER)")
set.labels <- rev(set.labels)

#Prep plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/GenePositiveControlSuppFig/",
          "GenePositiveControlSuppFig.dotplots.pdf",sep=""),
    height=8,width=12)
par(mar=c(0.5,21,3,1.5),bty="n",lend=1)
plot(x=c(-0.1,fold.max+0.1),y=c(0,nrow(x)),type="n",
     xaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")

#Gridlines & X-axis
abline(v=seq(0,fold.max,0.5),col="gray90",lwd=0.7)
abline(v=0:fold.max,col="gray85")
abline(v=0)
abline(v=1,lwd=2)
axis(3,at=seq(0,fold.max,0.5),col="gray50",tck=-0.01,labels=NA)
axis(3,at=0:fold.max,tck=-0.02,labels=NA)
axis(3,at=0:fold.max,line=-0.4,tick=F)
mtext(3,line=1.7,text="Fold-Enrichment versus Random Chance")

#Iterate and plot stats
sapply(1:nrow(x),function(i){
  plotVals <- as.numeric(x[i,-1])
  if(x[i,]$set!="SKIP"){
    #Account for NA fold-enrichments
    if(!is.na(plotVals[10]) & !is.na(plotVals[11])){
      #Account for fold-enrichments > fold.max
      if(plotVals[10]>fold.max){
        plotVals[10] <- fold.max
        plotVals[12] <- NA
        plotVals[13] <- NA
      }
      if(plotVals[11]>fold.max){
        plotVals[11] <- fold.max
        plotVals[14] <- NA
        plotVals[15] <- NA
      }
      #Shading
      rect(xleft=par("usr")[1],xright=par("usr")[2],
           ybottom=i-0.85,ytop=i-0.15,
           border=NA,col=adjustcolor("gray60",alpha=0.15))
      segments(x0=c(1,1),x1=plotVals[c(10,11)],
               y0=c(i-0.35,i-0.65),y1=c(i-0.35,i-0.65),
               lwd=6,col=adjustcolor(c("red","blue"),alpha=0.3))
      #CIs
      segments(x0=plotVals[c(12,14)],x1=plotVals[c(13,15)],
               y0=c(i-0.35,i-0.65),y1=c(i-0.35,i-0.65),
               lwd=1,col="gray20")
      #Point estimates
      points(x=plotVals[c(10,11)],y=c(i-0.35,i-0.65),
             pch=21,bg=c("red","blue"))
      #Significance marks
      if(plotVals[16]<0.001){
        axis(4,at=i-0.35,line=-0.9,tick=F,las=2,labels="***",col.axis="red")
      }else{
        if(plotVals[16]<0.01){
          axis(4,at=i-0.35,line=-0.9,tick=F,las=2,labels="**",col.axis="red")
        }else{
          if(plotVals[16]<0.05){
            axis(4,at=i-0.35,line=-0.9,tick=F,las=2,labels="*",col.axis="red")
          }
        }
      }
      if(plotVals[17]<0.001){
        axis(4,at=i-0.65,line=-0.9,tick=F,las=2,labels="***",col.axis="blue")
      }else{
        if(plotVals[17]<0.01){
          axis(4,at=i-0.65,line=-0.9,tick=F,las=2,labels="**",col.axis="blue")
        }else{
          if(plotVals[17]<0.05){
            axis(4,at=i-0.65,line=-0.9,tick=F,las=2,labels="*",col.axis="blue")
          }
        }
      }
      #Rounded fold-changes
      if((plotVals[10]>=1 & plotVals[10]<=(fold.max+1)/2) |
         (plotVals[10]>=0 & plotVals[10]<=0.5)){
        text(x=plotVals[10],y=i-0.15,pos=4,font=2,cex=0.8,
             col="red",labels=format(round(as.numeric(x[i,11]),2),nsmall=2))
      }else{
        text(x=plotVals[10],y=i-0.15,pos=2,font=2,cex=0.8,
             col="red",labels=format(round(as.numeric(x[i,11]),2),nsmall=2))
      }
      if((plotVals[11]>=1 & plotVals[11]<=(fold.max+1)/2) |
         (plotVals[11]>=0 & plotVals[11]<=0.5)){
        text(x=plotVals[11],y=i-0.85,pos=4,font=2,cex=0.8,
             col="blue",labels=format(round(as.numeric(x[i,12]),2),nsmall=2))
      }else{
        text(x=plotVals[11],y=i-0.85,pos=2,font=2,cex=0.8,
             col="blue",labels=format(round(as.numeric(x[i,12]),2),nsmall=2))
      }
    }else{
      #Shading
      rect(xleft=par("usr")[1],xright=par("usr")[2],
           ybottom=i-0.85,ytop=i-0.15,
           border=NA,col=adjustcolor("gray60",alpha=0.15))
      #Shadow point estimates
      points(x=c(1,1),y=c(i-0.35,i-0.65),
             pch=21,col="gray40",bg="gray80",
             lty=2,lwd=2)
      #NA labels
      axis(4,at=c(i-0.3,i-0.7),line=-0.9,tick=F,las=2,cex.axis=0.7,
           labels=c("NA","NA"),col.axis="gray50")
    }

    #Plot set label
    if(is.na(plotVals[2])){
      plotVals[2] <- 0
    }
    axis(2,at=i-0.5,tick=F,line=-0.9,las=2,cex.axis=0.8,
         labels=paste(set.labels[i]," [n=",
                      prettyNum(plotVals[2],big.mark=","),
                      " genes]",sep=""))
  }
})

#Close device
dev.off()


