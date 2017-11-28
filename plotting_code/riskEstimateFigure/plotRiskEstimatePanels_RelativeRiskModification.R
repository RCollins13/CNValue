#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to plot panels for risk estimate figure

#####Set parameters & load libraries
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.LOCI <- c("#1A5B22","#FFCB00","#34BFAD")
nsamps <- c(63629,57760,35693,22119,13361)
names(nsamps) <- c("GERM","NEURO","NDD","PSYCH","NONN")

#####PLOT A: incidence per pheno by category
#Read data
x <- lapply(list("largeSegment","gene","regBlock"),function(class){
  dat <- read.table(paste(WRKDIR,"plot_data/riskEstimatesFigure/",
                          class,".incidence.txt",sep=""),header=F)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  colnames(dat) <- c("DEL","DUP")
  return(dat)
})
names(x) <- c("largeSegment","gene","regBlock")

#Helper function to plot bars
plotBars <- function(df,ymax=NULL){
  #Set y-max
  if(is.null(ymax)){
    ymax <- 1.02*max(df)
  }

  #Prep plot area
  par(mar=c(1,2,1,0.5),bty="n")
  plot(x=c(0,(2*nrow(df))+2),y=c(0,ymax),type="n",
       yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Add gridlines
  abline(h=seq(0,floor(100*ymax)/100,0.01),col=cols.CTRL[3])

  #Plot bars
  rect(xleft=c((1:nrow(df))-1,seq(nrow(df)+1,(2*nrow(df)))+1),
       xright=c((1:nrow(df)),seq(nrow(df)+1,(2*nrow(df)))+2),
       ybottom=0,ytop=c(df[,1],df[,2]),
       col=c(cols.CTRL[1],cols.GERM[1],
             rep(cols.NEURO[1],3),cols.SOMA[1]),lwd=0.7)
  rect(xleft=c((1:nrow(df))-1,seq(nrow(df)+1,(2*nrow(df)))+1)[-c(1,7)],
       xright=c((1:nrow(df)),seq(nrow(df)+1,(2*nrow(df)))+2)[-c(1,7)],
       ybottom=0,ytop=c(rep(df[1,1],nrow(df)-1),
                        rep(df[1,2],nrow(df)-1)),
       col=cols.CTRL[2],lwd=0.7)

  #Plot axis
  axis(2,at=seq(0,floor(100*ymax)/100,0.02),labels=NA,lwd=0.7)
  axis(2,at=seq(0,floor(100*ymax)/100,0.02),line=-0.4,tick=F,las=2,
       labels=paste(seq(0,100*ymax,2),"%",sep=""))

  #Plot cleanup lines
  abline(h=0)
  rect(xleft=6.5,xright=7.5,
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border="white",col="white")
}

#Plot bars
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/riskEstimatesFigure/",
          "riskEstimateBarplots.pdf",sep=""),
    height=3,width=2)
par(mfrow=c(3,1))
plotBars(x$largeSegment,ymax=0.08)
plotBars(x$gene,ymax=0.08)
plotBars(x$regBlock,ymax=0.08)
dev.off()





#####Code to plot RRs & excess risk by locus & color by class
#Helper function to read data & return list of RRs
readRRs <- function(CNV){
  ###Read & clean segment data
  #Read table
  segs <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                           "1_2_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return RR & excess risk per signif pheno
  seg.vals <- t(apply(segs,1,function(vals){
    #Get segment ID
    ID <- paste(vals[1],vals[2],vals[3],sep="_")
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return RR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"p",sep="_"),colnames(segs))
      if(vals[idx]<0.05/12480){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(segs))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(segs))
        ncase.CNV <- vals[case.idx]
        nctrl.CNV <- vals[ctrl.idx]
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        ncase.noCNV <- ncase.total-ncase.CNV
        nctrl.noCNV <- nctrl.total-nctrl.CNV
        RR <- (ncase.CNV/(ncase.CNV+nctrl.CNV))/(ncase.noCNV/(ncase.noCNV+nctrl.noCNV))

        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(RR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant RRs
    if(!all(is.na(vals[,1]))){
      RR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      RR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(RR.final,risk.final,ID))
   }))
  seg.vals <- as.data.frame(seg.vals)
  colnames(seg.vals) <- c("RR","risk","ID")
  seg.vals$class <- "seg"
  seg.vals[,1:2] <- apply(seg.vals[,1:2],2,as.numeric)
  seg.vals <- seg.vals[which(!is.na(seg.vals$RR) & !is.infinite(seg.vals$RR)),]

  ###Read & clean gene data
  #Read table
  genes <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                           "3_4_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return RR & excess risk per signif pheno
  gene.vals <- t(apply(genes,1,function(vals){
    #Get genement ID
    ID <- vals[1]
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return RR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"q",sep="_"),colnames(genes))
      if(vals[idx]<0.05){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(genes))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(genes))
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        RR <- ((vals[case.idx]+1)/(vals[ctrl.idx]+(ncase.total/nctrl.total)))/
          ((ncase.total-vals[case.idx]+1E-6)/(nctrl.total-vals[ctrl.idx]+1E-6))
        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(RR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant RRs
    if(!all(is.na(vals[,1]))){
      RR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      RR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(RR.final,risk.final,ID))
  }))
  gene.vals <- as.data.frame(gene.vals)
  colnames(gene.vals) <- c("RR","risk","ID")
  gene.vals$class <- "gene"
  gene.vals[,1:2] <- apply(gene.vals[,1:2],2,as.numeric)
  gene.vals <- gene.vals[which(!is.na(gene.vals$RR) & !is.infinite(gene.vals$RR)),]

  ###Read & clean reg block data
  #Read table
  blocks <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                            "5_6_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return RR & excess risk per signif pheno
  block.vals <- t(apply(blocks,1,function(vals){
    #Get blockment ID
    ID <- paste(vals[1:3],collapse="_")
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return RR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"minimum_p",sep="_"),colnames(blocks))
      if(vals[idx]<0.05){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(blocks))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(blocks))
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        RR <- ((vals[case.idx]+1)/(vals[ctrl.idx]+(ncase.total/nctrl.total)))/
          ((ncase.total-vals[case.idx]+1E-6)/(nctrl.total-vals[ctrl.idx]+1E-6))
        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(RR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant RRs
    if(!all(is.na(vals[,1]))){
      RR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      RR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(RR.final,risk.final,ID))
  }))
  block.vals <- as.data.frame(block.vals)
  colnames(block.vals) <- c("RR","risk","ID")
  block.vals$class <- "reg"
  block.vals[,1:2] <- apply(block.vals[,1:2],2,as.numeric)
  block.vals <- block.vals[which(!is.na(block.vals$RR) & !is.infinite(block.vals$RR)),]

  ###Prepare final results & reorder
  dat.out <- rbind(seg.vals,gene.vals,block.vals)
  dat.out[,1:2] <- apply(dat.out[,1:2],2,as.numeric)
  dat.out <- dat.out[order(dat.out$RR,decreasing=T),]
  return(dat.out)
}

#Helper function to make ranked barplots of all RRs and risk by class
rankBarplot <- function(vals,classes,
                        ymax=NULL,xlab="Loci",
                        meanBar="black"){
  #Convert classes to colors
  cols <- classes
  cols <- gsub("seg",cols.LOCI[1],cols)
  cols <- gsub("gene",cols.LOCI[2],cols)
  cols <- gsub("reg",cols.LOCI[3],cols)

  #Set ymax
  if(is.null(ymax)){
    ymax <- 1.02*max(vals)
  }

  #Prep plot area
  par(mar=c(1.2,2.8,1,1),bty="n")
  plot(x=c(0,length(vals)+1),y=c(0,ymax),type="n",
       yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")

  #Add gridlines
  abline(h=seq(0,ceiling(par("usr")[2]),0.5),lwd=0.5,col=cols.CTRL[3])
  abline(h=seq(0,ceiling(par("usr")[2])),lwd=0.75,col=cols.CTRL[2])

  #Add axes & labels
  mtext(1,text=paste(xlab," (n=",length(vals),")",sep=""),line=0.1)
  axis(2,at=seq(0,ceiling(par("usr")[2]),0.5),tck=-0.01,col=cols.CTRL[1],labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[2])),labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[2])),tick=F,las=2,line=-0.3,cex.axis=0.8,
       labels=2^seq(0,ceiling(par("usr"))[2]))
  mtext(2,text=expression(paste(RR,phantom(x),(log[2]-scaled))),line=1.5)

  #Plot rectangles
  rect(xleft=seq(0,length(vals)-1),xright=seq(1,length(vals)),
       ybottom=0,ytop=vals,border=cols,col=cols)

  #Misc clean up
  abline(h=0)
  abline(h=mean(vals),col=meanBar)
  print(2^mean(vals))
}

#Helper function to make swarmplots of RRs

#Helper function to make scatterplots of RR vs excess risk



#Read data
DEL <- readRRs("DEL")
DUP <- readRRs("DUP")

#Plot ranked barplots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/riskEstimatesFigure/",
          "RRs_per_locus.barplots.pdf",sep=""),
    height=4.8,width=4)
par(mfrow=c(2,1))
rankBarplot(log2(DEL$RR),DEL$class,ymax=8,
            xlab="DEL-significant Loci",meanBar="red")
rankBarplot(log2(DUP$RR),DUP$class,ymax=8,
            xlab="DUP-significant Loci",meanBar="blue")
dev.off()


boxplot(log2(DEL[which(DEL$class=="seg"),1]),
        log2(DEL[which(DEL$class=="gene"),1]),
        log2(DEL[which(DEL$class=="reg"),1]),
        col=cols.LOCI,ylim=c(0,8))

boxplot(log2(DUP[which(DEL$class=="seg"),1]),
        log2(DUP[which(DEL$class=="gene"),1]),
        log2(DUP[which(DEL$class=="reg"),1]),
        col=cols.LOCI,ylim=c(0,8))

boxplot(DEL[which(DEL$class=="seg"),2],
        DEL[which(DEL$class=="gene"),2],
        DEL[which(DEL$class=="reg"),2],
        col=cols.LOCI)

boxplot(DUP[which(DEL$class=="seg"),2],
        DUP[which(DEL$class=="gene"),2],
        DUP[which(DEL$class=="reg"),2],
        col=cols.LOCI)





