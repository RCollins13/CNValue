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
require(beeswarm)
require(MASS)

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





#####Code to plot ORs & excess risk by locus & color by class
#Helper function to read data & return list of ORs
readORs <- function(CNV){
  ###Read & clean segment data
  #Read table
  segs <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                           "1_2_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return OR & excess risk per signif pheno
  seg.vals <- t(apply(segs,1,function(vals){
    #Get segment ID
    ID <- paste(vals[1],vals[2],vals[3],sep="_")
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return OR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"p",sep="_"),colnames(segs))
      if(vals[idx]<0.05/12480){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(segs))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(segs))
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        OR <- ((vals[case.idx]+0.5)/(vals[ctrl.idx]+0.5))/
          ((ncase.total-vals[case.idx]+0.5)/(nctrl.total-vals[ctrl.idx]+0.5))
        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(OR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant ORs
    if(!all(is.na(vals[,1]))){
      OR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      OR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(OR.final,risk.final,ID))
  }))
  seg.vals <- as.data.frame(seg.vals)
  colnames(seg.vals) <- c("OR","risk","ID")
  seg.vals$class <- "seg"
  seg.vals[,1:2] <- apply(seg.vals[,1:2],2,as.numeric)
  seg.vals <- seg.vals[which(!is.na(seg.vals$OR) & !is.infinite(seg.vals$OR)),]

  ###Read & clean gene data
  #Read table
  genes <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                            "3_4_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return OR & excess risk per signif pheno
  gene.vals <- t(apply(genes,1,function(vals){
    #Get genement ID
    ID <- vals[1]
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return OR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"q",sep="_"),colnames(genes))
      if(vals[idx]<0.05){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(genes))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(genes))
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        OR <- ((vals[case.idx]+0.5)/(vals[ctrl.idx]+0.5))/
          ((ncase.total-vals[case.idx]+0.5)/(nctrl.total-vals[ctrl.idx]+0.5))
        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(OR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant ORs
    if(!all(is.na(vals[,1]))){
      OR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      OR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(OR.final,risk.final,ID))
  }))
  gene.vals <- as.data.frame(gene.vals)
  colnames(gene.vals) <- c("OR","risk","ID")
  gene.vals$class <- "gene"
  gene.vals[,1:2] <- apply(gene.vals[,1:2],2,as.numeric)
  gene.vals <- gene.vals[which(!is.na(gene.vals$OR) & !is.infinite(gene.vals$OR)),]

  ###Read & clean reg block data
  #Read table
  blocks <- read.table(paste(WRKDIR,"plot_data/suppTables/suppTables_",
                             "5_6_",CNV,".txt",sep=""),header=T)
  #Iterate over sites & return OR & excess risk per signif pheno
  block.vals <- t(apply(blocks,1,function(vals){
    #Get blockment ID
    ID <- paste(vals[1:3],collapse="_")
    #Convert vals to numeric
    vals <- as.numeric(vals)
    #Iterate over phenotypes & return OR & excess risk for significant phenos
    vals <- t(sapply(c("GERM","NEURO","NDD","PSYCH","NONN"),function(pheno){
      idx <- grep(paste(pheno,CNV,"minimum_p",sep="_"),colnames(blocks))
      if(vals[idx]<0.05){
        case.idx <- grep(paste(pheno,CNV,"count",sep="_"),colnames(blocks))
        ctrl.idx <- grep(paste("CTRL",CNV,"count",sep="_"),colnames(blocks))
        ncase.total <- nsamps[which(names(nsamps)==pheno)]
        nctrl.total <- 38628
        OR <- ((vals[case.idx]+0.5)/(vals[ctrl.idx]+0.5))/
          ((ncase.total-vals[case.idx]+0.5)/(nctrl.total-vals[ctrl.idx]+0.5))
        risk <- (vals[case.idx]/ncase.total)-(vals[ctrl.idx]/nctrl.total)
        return(c(OR,risk))
      }else{
        return(c(NA,NA))
      }
    }))
    #Weighted average of significant ORs
    if(!all(is.na(vals[,1]))){
      OR.final <- sum(vals[,1]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,1]))])
    }else{
      OR.final <- NA
    }
    #Weighted average of significant risk estimates
    if(!all(is.na(vals[,2]))){
      risk.final <- sum(vals[,2]*nsamps,na.rm=T)/sum(nsamps[which(!is.na(vals[,2]))])
    }else{
      risk.final <- NA
    }

    #Return values
    return(c(OR.final,risk.final,ID))
  }))
  block.vals <- as.data.frame(block.vals)
  colnames(block.vals) <- c("OR","risk","ID")
  block.vals$class <- "reg"
  block.vals[,1:2] <- apply(block.vals[,1:2],2,as.numeric)
  block.vals <- block.vals[which(!is.na(block.vals$OR) & !is.infinite(block.vals$OR)),]

  ###Prepare final results & reorder
  dat.out <- rbind(seg.vals,gene.vals,block.vals)
  dat.out[,1:2] <- apply(dat.out[,1:2],2,as.numeric)
  dat.out <- dat.out[order(dat.out$OR,decreasing=T),]
  return(dat.out)
}

#Helper function to make ranked barplots of all ORs and risk by class
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
  plot(x=c(0,length(vals)+1),y=c(0,1.02*ymax),type="n",
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
  mtext(2,text=expression(paste(OR,phantom(x),(log[2]-scaled))),line=1.5)

  #Plot rectangles
  rect(xleft=seq(0,length(vals)-1),xright=seq(1,length(vals)),
       ybottom=0,ytop=vals,border=NA,col=cols)

  #Add asterixes above plot for sites that extend above ymax
  sapply(which(vals>ymax),function(i){
    axis(3,at=i,line=-1.5,tick=F,labels="*",col.axis=cols[i])
  })

  #Misc clean up
  abline(h=0)
  abline(h=mean(vals),col=meanBar)
  print(2^mean(vals))
}

#Helper function to make swarmplots of ORs
swarmORs <- function(vals,classes,ymax=NULL){
  #Convert classes to colors
  cols <- classes
  cols <- gsub("seg",cols.LOCI[1],cols)
  cols <- gsub("gene",cols.LOCI[2],cols)
  cols <- gsub("reg",cols.LOCI[3],cols)

  #Log2-transform vals
  vals <- log2(vals)

  #Set ymax
  if(is.null(ymax)){
    ymax <- 1.02*max(vals)
  }

  #Round values down to ymax
  vals[which(vals>ymax)] <- ymax

  #Prep plot area
  par(mar=c(1.2,2.8,1,1),bty="n")
  plot(x=c(0,3),y=c(0,1.02*ymax),type="n",
       yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")

  #Add gridlines
  abline(h=seq(0,ceiling(par("usr")[4]),0.5),lwd=0.5,col=cols.CTRL[3])
  abline(h=seq(0,ceiling(par("usr")[4])),lwd=0.75,col=cols.CTRL[2])

  #Add axes & labels
  # sapply(1:3,function(i){
  #   axis(1,at=i-0.5,tick=F,line=0.2,
  #        labels=paste(c("Seg.","Gene","Reg.")[i],"\n(n=",
  #                     prettyNum(length(which(classes==c("seg","gene","reg")[i])),big.mark=","),
  #                     ")",sep=""))
  # })
  axis(2,at=seq(0,ceiling(par("usr")[4]),0.5),tck=-0.01,col=cols.CTRL[1],labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[4])),labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[4])),tick=F,las=2,line=-0.3,cex.axis=0.8,
       labels=2^seq(0,ceiling(par("usr"))[4]))
  # mtext(2,text=expression(paste(OR,phantom(x),(log[2]-scaled))),line=1.5)

  #Swarmplots
  sapply(1:3,function(i){
    vals.c <- vals[which(classes==c("seg","gene","reg")[i])]
    #Swarmplot
    beeswarm(vals.c,add=T,xlab="",ylab="",xaxt="n",yaxt="n",
             at=i-0.5,corral="random",corralWidth=0.8,
             pch=19,col=cols.LOCI[i],cex=0.5)
  })

  # #Reference marker line
  # means <- sapply(c("seg","gene","reg"),function(class){
  #   return(mean(vals[which(classes==class)]))
  # })
  # segments(x0=0.5:1.5,x1=1.5:2.5,
  #          y0=means[1:2],y1=means[2:3],
  #          lwd=6,col=cols.CTRL[3])

  #Mean & IQR
  sapply(1:3,function(i){
    vals.c <- vals[which(classes==c("seg","gene","reg")[i])]
    #IQR
    quants <- quantile(vals.c,p=c(0.25,0.75))
    segments(x0=i-0.65,x1=i-0.35,y0=quants,y1=quants,lwd=2)
    # segments(x0=i-0.5,x1=i-0.5,y0=quants[1],y1=quants[2],lwd=2)
    #Mean
    points(x=i-0.5,y=mean(vals.c),pch=23,col="black",bg="white",lwd=2,cex=1.5)
  })
}

#Helper function to make scatterplots of OR vs excess risk
RiskVsOR <- function(OR,risk,classes,xmax=NULL,ymax=NULL,xstep=0.00025){
  #Convert classes to colors
  cols <- classes
  cols <- gsub("seg",cols.LOCI[1],cols)
  cols <- gsub("gene",cols.LOCI[2],cols)
  cols <- gsub("reg",cols.LOCI[3],cols)

  #Log2-transform ORs
  OR <- log2(OR)

  #Set xmax and ymax
  if(is.null(ymax)){
    ymax <- 1.02*max(OR)
  }
  if(is.null(xmax)){
    xmax <- 1.02*max(risk)
  }

  #Round values down to max
  OR[which(OR>ymax)] <- ymax
  risk[which(risk>xmax)] <- xmax

  #Prep plot area
  par(mar=c(2.8,2.8,1,1),bty="n")
  plot(x=c(0,1.02*xmax),y=c(0,1.02*ymax),type="n",
       xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")

  #Add gridlines
  abline(h=seq(0,ceiling(par("usr")[4]),0.5),lwd=0.5,col=cols.CTRL[3])
  abline(h=seq(0,ceiling(par("usr")[4])),lwd=0.75,col=cols.CTRL[2])
  abline(v=seq(0,par("usr")[2]+0.001,xstep),lwd=0.5,col=cols.CTRL[3])
  abline(v=seq(0,par("usr")[2]+0.001,2*xstep),lwd=0.75,col=cols.CTRL[2])

  #Add axes & labels
  # sapply(1:3,function(i){
  #   axis(1,at=i-0.5,tick=F,line=0.2,
  #        labels=paste(c("Seg.","Gene","Reg.")[i],"\n(n=",
  #                     prettyNum(length(which(classes==c("seg","gene","reg")[i])),big.mark=","),
  #                     ")",sep=""))
  # })
  axis(2,at=seq(0,ceiling(par("usr")[4]),0.5),tck=-0.01,col=cols.CTRL[1],labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[4])),labels=NA)
  axis(2,at=seq(0,ceiling(par("usr")[4])),tick=F,las=2,line=-0.3,cex.axis=0.8,
       labels=2^seq(0,ceiling(par("usr"))[4]))
  mtext(2,text=expression(paste(OR,phantom(x),(log[2]-scaled))),line=1.5)
  axis(1,at=seq(0,par("usr")[2]+0.001,xstep),tck=-0.01,col=cols.CTRL[1],labels=NA)
  axis(1,at=seq(0,par("usr")[2]+0.001,2*xstep),labels=NA)
  axis(1,at=seq(0,par("usr")[2]+0.001,2*xstep),tick=F,line=-0.4,cex.axis=0.8,
       labels=paste(100*seq(0,par("usr")[2]+0.001,2*xstep),"%",sep=""))
  mtext(1,text="Excess case risk (% cases - % controls)",line=1.5)

  #Convex hulls
  sapply(c("seg","gene","reg"),function(class){
    col.class <- cols.LOCI[which(c("seg","gene","reg")==class)]
    k <- chull(risk[which(classes==class)],
          OR[which(classes==class)])
    polygon(x=risk[which(classes==class)][k],
            y=OR[which(classes==class)][k],
            border=adjustcolor(col.class,alpha=0.7),lty=2,
            col=adjustcolor(col.class,alpha=0.3))
  })

  #Scatterplot
  points(x=risk,y=OR,pch=19,col=cols)

  #Means
  sapply(c("seg","gene","reg"),function(class){
    col.class <- cols.LOCI[which(c("seg","gene","reg")==class)]
    m.x <- mean(risk[which(classes==class)])
    m.y <- mean(OR[which(classes==class)])
    points(x=m.x,y=m.y,pch=21,cex=1.5,col="black",bg="white")
  })
}



#Read data
DEL <- readORs("DEL")
DUP <- readORs("DUP")

#Plot ranked barplots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/riskEstimatesFigure/",
          "ORs_per_locus.barplots.pdf",sep=""),
    height=4.8,width=3)
par(mfrow=c(2,1))
rankBarplot(log2(DEL$OR),DEL$class,ymax=7,
            xlab="DEL Loci",meanBar="red")
rankBarplot(log2(DUP$OR),DUP$class,ymax=7,
            xlab="DUP Loci",meanBar="blue")
dev.off()

#Plot swarmplots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/riskEstimatesFigure/",
          "ORs_per_locus.swarmplots.pdf",sep=""),
    height=4.8,width=1.8)
par(mfrow=c(2,1))
swarmORs(DEL$OR,DEL$class,ymax=7)
swarmORs(DUP$OR,DUP$class,ymax=7)
dev.off()


boxplot(log2(DEL[which(DEL$class=="seg"),1]),
        log2(DEL[which(DEL$class=="gene"),1]),
        log2(DEL[which(DEL$class=="reg"),1]),
        col=cols.LOCI,ylim=c(0,8))

boxplot(log2(DUP[which(DUP$class=="seg"),1]),
        log2(DUP[which(DUP$class=="gene"),1]),
        log2(DUP[which(DUP$class=="reg"),1]),
        col=cols.LOCI,ylim=c(0,8))

boxplot(DEL[which(DEL$class=="seg"),2],
        DEL[which(DEL$class=="gene"),2],
        DEL[which(DEL$class=="reg"),2],
        col=cols.LOCI)

boxplot(DUP[which(DEL$class=="seg"),2],
        DUP[which(DEL$class=="gene"),2],
        DUP[which(DEL$class=="reg"),2],
        col=cols.LOCI)





