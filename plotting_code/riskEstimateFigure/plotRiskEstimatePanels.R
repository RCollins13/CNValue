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
    # points(x=i-0.5,y=mean(vals.c),pch=23,col="black",bg="white",lwd=2,cex=1.5)
    points(x=i-0.5,y=mean(vals.c),pch=23,cex=1.5,col="black",bg="white")
    points(x=i-0.5,y=mean(vals.c),pch=23,cex=0.75,lwd=0.7,col="black",bg=cols.LOCI[i])
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
  par(mar=c(2.8,2.8,1.5,1.5),bty="n")
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
  mtext(1,text="Excess case risk",line=1.5)

  #Scatterplot
  points(x=risk,y=OR,pch=19,cex=0.5,col=cols)

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

  #Means
  sapply(c("seg","gene","reg"),function(class){
    col.class <- cols.LOCI[which(c("seg","gene","reg")==class)]
    m.x <- mean(risk[which(classes==class)])
    m.y <- mean(OR[which(classes==class)])
    points(x=m.x,y=m.y,pch=23,cex=1.5,col="black",bg="white")
    points(x=m.x,y=m.y,pch=23,cex=0.75,lwd=0.7,col="black",bg=col.class)
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

#Plot scatterplots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/riskEstimatesFigure/",
          "ORs_vs_risk_per_locus.scatterplots.pdf",sep=""),
    height=5,width=2.5)
par(mfrow=c(2,1))
RiskVsOR(DEL$OR,DEL$risk,DEL$class,ymax=7,xmax=0.003,xstep=0.0005)
RiskVsOR(DUP$OR,DUP$risk,DUP$class,ymax=7,xmax=0.003,xstep=0.0005)
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


#####Test plots when prototyping the SSC proband-sibling filter
keep_DEL_genes_bonf <- c("AC132872.2","ADAMTS18","ADRBK2","AGAP1","AHCTF1","ASTN2","ATP9B",
                    "AUTS2","BACH2","BNC2","C14orf132","C7orf10","CADM1","CALN1","CD82",
                    "CDH18","CDH19","CDKAL1","CHL1","CNBD1","CSGALNACT1","CSMD1","CSMD3",
                    "CSNK1D","CTNND2","DLG2","DLGAP1","DLGAP2","DOCK8","DTWD2","EHMT1",
                    "EML6","ERBB4","ERICH1","FAM46C","FARS2","FBN2","FBXL4","FGD4","FRG1",
                    "GBE1","GLIS3","GOLPH3","GPC5","GRIP1","GTDC1","HS3ST2","IMMP2L","ISPD",
                    "KANK1","KCNV2","KIAA0020","LARGE","LPIN2","LPP","LTBP1","MAD1L1","MBD5",
                    "MED13L","MEF2C","MEI4","MME","MNAT1","NAALADL2","NEGR1","NKAIN2","NRXN1",
                    "ORC4","OVGP1","PACRG","PCDH9","PDE4D","PDZD2","PLD5","PPP2R2B","PTPRD",
                    "PTPRT","RAB3GAP2","RAPGEF2","RXRA","SDK1","SGCZ","SIX2","SLC16A3","SLC1A1",
                    "SLC2A1","SMYD3","SOX5","SPAG16","STK17B","SUMF1","SUPT3H","SYT1","TBC1D22A",
                    "TBL1XR1","TCF4","THADA","TRANK1","TRIM32","VIPR2","WNT11","WSCD1","ZEB2",
                    "ZFP42","ZNF385D","ZNF804B")
keep_DEL_genes_pheno <- c("AC132872.2","ASTN2","AUTS2","B3GNTL1","BACH2","CD82","CDC42","CDKAL1",
                          "CSNK1D","CTNND2","DLG2","ERICH1","FAM19A5","FAM46C","FARS2","GBE1",
                          "GNAT3","GOLPH3","GTDC1","HS3ST2","ISPD","KCNV2","KIAA0020","LEF1",
                          "LRRFIP2","LYRM4","MBD5","MED13L","MEF2C","NRXN1","ORC4","PDZD2","RAB3GAP2",
                          "RFX3","SAMSN1","SIX2","SLC16A3","SLC1A1","SLC2A1","SLCO1B1","SMARCA2","SOX5",
                          "STK17B","TBL1XR1","TCF4","TPD52L3","TRIM32","VIPR2","WNT11","ZEB2","ZMYND11")
keep_DEL_genes <- intersect(keep_DEL_genes_pheno,keep_DEL_genes_bonf)
DEL <- DEL[which(DEL$ID %in% keep_DEL_genes),]
gt_DEL_genes <- c("NRXN1","MBD5","ORC4","ARHGEF10","GBE1","ISPD","MALL","ATP9B",
                  "CDC42","COL6A2","DLGAP2","ERBB4","FAM196A","FAM19A5","FBXO25",
                  "GLIS3","KCNV2","KIAA0020","PSMD7","SLCO1B1","SOX5","TNRC6B")
lt_DEL_genes <- c("CDH18","DLG2","DOCK8","DTWD2","ZNF14")
keep_DUP_genes_bonf <- c("A1CF","ABCC1","ABCC6","AC092850.1","ACSM4","ACTR3B","ACTRT2","ADAMTS17",
                    "ADAMTS2","ADCYAP1","AHI1","AHR","AJAP1","AL450307.1","ALG10B","ANGPT1",
                    "ANKFY1","ANKRD30B","AP003062.1","APBA2","ASAH2B","ATAD5","ATG5","ATP2C2",
                    "AUTS2","BPTF","C16orf62","C16orf72","C17orf51","C17orf58","C18orf42",
                    "C2orf27A","C5orf42","CACNA1B","CACNA2D1","CALN1","CCDC50","CEP85L",
                    "CHL1","CLIC5","CNTN4","CNTN5","CNTNAP4","CRLF3","CSMD1","CSMD3","CTNNA2",
                    "CTNND2","CWH43","CYB5D2","DEFB115","DIP2C","DLG1","DMRT1","DMRT2","DOK6",
                    "DPP6","EMB","ERC1","ERC2","FAM110C","FAM150A","FAM194B","FAM19A5","FBRSL1",
                    "FBXL7","FOXD2","FRG2C","GALNT11","GALNTL6","GHR","GJA5","GK5","GLDC","GLRA1",
                    "GMDS","GPBP1","GRM3","HERC2","IMMP2L","KANK1","KDM4C","KIF26B","KIRREL3",
                    "KMT2C","LAMA1","LINC00955","LINGO1","LTBP1","LYZL1","MACROD2","MALL","MARCH1",
                    "MEI4","METTL4","MPRIP","MRGPRE","MRPL1","MSR1","MTRNR2L1","MTX2","MUC8",
                    "NCALD","NEBL","NEK10","NINL","NPHP1","NRXN1","NTN1","NUDT12","NUP155","OR4A5",
                    "PAPPA","PARD3","PCNX","PDE5A","PFKP","PLEKHM1","PLN","PODXL","PRIM2","PTPRJ",
                    "PTPRM","PTPRN2","PVRL1","RABL2A","RAF1","RASA3","RBFOX1","RHOU","RREB1","RRM1",
                    "RYR2","SAMD5","SCAPER","SDK1","SEMA5A","SGCZ","SHANK3","SHCBP1","SLIT3","SLITRK6",
                    "SMOC2","SMYD3","SNTG1","SPOCK3","SRBD1","ST18","ST3GAL3","SUCLG1","SYK","TAS2R1",
                    "TBCE","TCERG1L","TFRC","TGM6","THBS2","TMEM179","TMEM18","TPO","TRIM48","TRIML1",
                    "TRIML2","TTC40","TUSC1","UBAP1","UBBP4","USP7","UTS2B","VIPR2","VWDE","XPO6","YES1",
                    "ZFAND3","ZFP42","ZMYND11","ZNF195","ZNF385D","ZNF717")
keep_DUP_genes_pheno <- c("ABCC6","AC093802.1","ACP1","ACTRT2","ADAMTS17","ADAMTS2","AHI1","AL450307.1",
                          "ANKFY1","ANKRD30A","ANKRD30B","AP003062.1","APBA2","ARHGAP28","ATAD5","ATG5",
                          "AUTS2","BICC1","BPTF","C11orf44","C14orf144","C16orf72","C17orf58","C20orf187",
                          "C5orf42","CACNA1B","CD8A","CLIC5","CLUH","CNGA1","COL23A1","COLEC12","CRLF3",
                          "CWH43","CYB5D2","DEFB115","DIP2C","DLG5","DMRT1","DMRT2","DOK6","EFCAB2","EMB",
                          "ESYT2","EXOC2","FAM110C","FAM150B","FAM19A5","FGF12","FRMPD2","GALNT11","GEM",
                          "GFRAL","GMDS","GRM3","HERC2","HUS1B","KIF26B","LINC00955","LRRC30","MACROD2",
                          "MALL","MARCH1","MEI4","MRGPRE","MSH3","MTRNR2L1","MTX2","MYLK3","NAMPTL","NBPF3",
                          "NCAPG2","NINL","NIPAL1","NPHP1","NTN1","NUDT12","ORC6","PARD3","PFKP","PLEKHM1",
                          "PPP2R3C","PTPRJ","RAB11FIP4","RAF1","RAP1GAP2","RASA3","RMND5A","RPSAP58","RXFP2",
                          "SAMD5","SGCZ","SH3YL1","SHANK3","SHCBP1","SMARCA2","SMOC2","SMYD3","SNCAIP","SNTG1",
                          "SPOCK3","SUCLG1","TAS2R1","TMEM179","TPO","TRIML1","TRIML2","TUSC1","TXNRD1",
                          "UBAP1","UNCX","USP7","VIPR2","VPS35","WDR60","YES1","ZFP42","ZNF14","ZNF195","ZNF675","ZNF681")
keep_DUP_genes <- intersect(keep_DUP_genes_pheno,keep_DUP_genes_bonf)
DUP <- DUP[which(DUP$ID %in% keep_DUP_genes),]
gt_DUP_genes <- c("HERC2","DMRT1","GJA5","NPHP1","ATAD5","CRLF3","NUP155","SGCZ",
                  "C5orf42","CTDP1","DEFB115","DLG1","DMRT2","DNAH5","PROP1","RASA3",
                  "SRBD1","USP7","ADAMTS17","ADAMTS2","AGRN","AL450307.1","APBA2",
                  "ATP2C2","C16orf45","C16orf72","CD8A","CHMP4B","CHST9","CNTNAP4",
                  "COLEC12","DOK6","FAM189A1","FAM53A","FRMPD2","GMDS","KIAA0020",
                  "LINC00955","MALL","NINL","RMND5A","RNF223","ROCK1","SMARCA2",
                  "SPAG16","SPOCK3","TMEM179","TMEM18","TPTE2","ZNF732")
lt_DUP_genes <- c("ADCYAP1","ANKRD30B","ARHGAP28","C18orf42","C1orf222","CNGA1",
                  "EFCAB2","EXOC2","GALNTL6","GLRA1","HUS1B","LRRC30","MRGPRE",
                  "NIPAL1","OR2V1","PARD3","PFKP","PPP2R3C","PRR26","RPSAP58","SHANK3",
                  "TRAT1","WDR60","YES1","ZNF195","ZNF675","ZNF681","ANKFY1","CYB5D2",
                  "EMB","GFRAL","OR4A15","OR4A16","TMEM40","TRIM48","TRIML1","TRIML2",
                  "VIPR2","ZMYND11","ZNF66","CACNA1B","RAF1","MTRNR2L1","SMYD3")
plot(x=c(0,6),y=c(0,7),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
axis(2,at=seq(0,ceiling(par("usr")[4]),0.5),tck=-0.01,col=cols.CTRL[1],labels=NA)
axis(2,at=seq(0,ceiling(par("usr")[4])),labels=NA)
axis(2,at=seq(0,ceiling(par("usr")[4])),tick=F,las=2,line=-0.3,cex.axis=0.8,
     labels=2^seq(0,ceiling(par("usr"))[4]))
mtext(2,text=expression(paste(OR,phantom(x),(log[2]-scaled))),line=1.5)
axis(1,at=0.5:5.5,labels=c("DEL\nPro>Sib","DEL\nPro=Sib\nOr no CNVs","DEL\nPro<Sib",
                           "DUP\nPro>Sib","DUP\nPro=Sib\nOr no CNVs","DUP\nPro<Sib"),
     tick=F,line=2)
#Mini func
plotVals <- function(vals,col,at){
  beeswarm(vals,col=col,add=T,at=at,pch=19,corral="random",corralWidth=0.8,xaxt="n",yaxt="n")
  boxplot(vals,add=T,at=at,col=NA,lwd=3,xaxt="n",yaxt="n")
}
#DEL gt
plotVals(log2(DEL[which(DEL$ID %in% gt_DEL_genes),1]),
         "red",0.5)
#DEL eq/NA
plotVals(log2(DEL[which(!(DEL$ID %in% gt_DEL_genes) & !(DEL$ID %in% lt_DEL_genes)),1]),
         "red",1.5)
#DEL lt
plotVals(log2(DEL[which(DEL$ID %in% lt_DEL_genes),1]),
         "red",2.5)
#DUP gt
plotVals(log2(DUP[which(DUP$ID %in% gt_DUP_genes),1]),
         "blue",3.5)
#DUP eq/NA
plotVals(log2(DUP[which(!(DUP$ID %in% gt_DUP_genes) & !(DUP$ID %in% lt_DUP_genes)),1]),
         "blue",4.5)
#DUP lt
plotVals(log2(DUP[which(DUP$ID %in% lt_DUP_genes),1]),
         "blue",5.5)
title("Comparison of Gene ORs by SSC CNV Evidence (Analysis-wide Bonferroni-sig AND just GERM/NEURO/NDD)")



