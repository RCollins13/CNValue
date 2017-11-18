#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate mirrored manhattan plots of major phenotype groups

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

###################################################
#####Helper function to read and transform p-values
###################################################
readData <- function(pheno,VF,filt){
  df <- read.table(paste(WRKDIR,"plot_data/figure2/",pheno,".",
                         VF,".",filt,".manhattan_pvals.bed",sep=""),header=F)
  names(df) <- c("chr","start","end","pDEL","pDUP")
  df$pDEL <- -log10(df$pDEL)
  df$pDUP <- -log10(df$pDUP)
  df$pos <- apply(df[,2:3],1,mean)
  return(df[,c(1,6,4:5)])
}

################################################
######Helper function to plot mirrored manhattan
################################################
mirrorManhattan <- function (df,adjusted=1E-8,ymax=NULL,col.even,col.odd,yaxis=T){
  #Gather list of unique contigs
  contigs <- unique(df[,1])
  contigs <- contigs[which(!(is.na(contigs)))]

  #Create index data frame
  indexes <- as.data.frame(t(sapply(contigs,function(chr){
    return(c(chr,0,max(df[which(df[,1]==chr),2])))
  })))
  indexes$sum <- cumsum(indexes[,3])
  indexes$bg <- rep(c(col.even,col.odd),ceiling(nrow(indexes)/2))[1:nrow(indexes)]

  #Create new plotting df with modified coordinates & color assignments
  df.plot <- as.data.frame(t(apply(df,1,function(row){
    #Determine del/dup colors
    if(as.numeric(row[3]>-log10(adjusted))){
      delCol <- "red"
    }else{
      delCol <- indexes[as.numeric(row[1]),5]
    }
    if(as.numeric(row[4]>-log10(adjusted))){
      dupCol <- "blue"
    }else{
      dupCol <- indexes[as.numeric(row[1]),5]
    }
    #Return vector
    return(c(row[1],
             as.numeric(row[2]+indexes[as.numeric(row[1]),4]-indexes[as.numeric(row[1]),3]),
             as.numeric(row[3]),
             as.numeric(row[4]),
             delCol,dupCol))
  })))
  df.plot[,2] <- as.numeric(as.character(df.plot[,2]))
  df.plot[,3] <- as.numeric(as.character(df.plot[,3]))
  df.plot[,4] <- as.numeric(as.character(df.plot[,4]))

  #Round values to ymax (if optioned)
  if(is.null(ymax)){
    ymax <- max(df.plot[,3:4],na.rm=T)
  }else{
    df.plot[which(df.plot[,3]>ymax),3] <- ymax
    df.plot[which(df.plot[,4]>ymax),4] <- ymax
  }

  #Set height of contig boxes
  boxht <- 0.05*ymax

  #Prepare plotting window
  par(mar=c(0.2,2,0.2,0.2))
  plot(x=c(0,max(indexes[,4])),y=c(-1.1*(ymax+boxht),1.1*(ymax+boxht)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Add background shading
  rect(xleft=rep(par("usr")[1],3),xright=rep(par("usr")[2],3),
       ybottom=c(par("usr")[3],log10(adjusted)-boxht,log10(0.05)-boxht),
       ytop=c(par("usr")[4],-log10(adjusted)+boxht,-log10(0.05)+boxht),
       border=NA,col=cols.CTRL[4:2])
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=-boxht,ytop=boxht,col="white")

  #Add y-axis
  if(yaxis==T){
    y.at <- seq(0,2*ceiling(par("usr")[4]),by=5)
    #Top y-axis
    axis(2,at=y.at+boxht,labels=NA)
    axis(2,at=y.at[-1]+boxht,tick=F,line=-0.3,labels=y.at[-1],las=2)
    #Bottom y-axis
    axis(2,at=-(y.at+boxht),labels=NA)
    axis(2,at=-(y.at[-1]+boxht),tick=F,line=-0.3,labels=y.at[-1],las=2)
  }

  #Add horizontal gridlines
  abline(v=indexes[,4],col="white",lwd=2)
  abline(h=c(-seq(0,ymax+10,5)-boxht,seq(0,ymax+10,5)+boxht),
         col=cols.CTRL[3],lwd=0.8)
  abline(h=c(log10(adjusted)-boxht,-log10(adjusted)+boxht),
         col=adjustcolor(c("red","blue"),alpha=0.5),lty=3)

  #Plots points
  points(df.plot[,2],df.plot[,4]+boxht,
         cex=0.5,pch=19,col=as.character(df.plot[,6]))
  points(df.plot[,2],-(df.plot[,3]+boxht),
         cex=0.5,pch=19,col=as.character(df.plot[,5]))

  #Adds chromosome labels
  sapply(1:length(contigs),function(i){
    #Rectangle
    rect(xleft=indexes[i,4]-indexes[i,3],xright=indexes[i,4],
         ybottom=-boxht,ytop=boxht,col="white")
    if(i<=13){
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.85)
    }else if(i<=18){
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.6,srt=90)
    }else{
      text(x=mean(c(indexes[i,4]-indexes[i,3],xright=indexes[i,4])),y=0,
           labels=contigs[i],font=4,cex=0.45,srt=90)
    }
  })

  # #Add dashed lines at ymax
  # abline(h=c(ymax+boxht,-(ymax+boxht)),lty=2)

  #Add cleanup box
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA,border="black",lwd=2)
}

########################################
#####Generate manhattan plots for figure
########################################
#Make plotting values list
plotVals <- list(c("GERM",cols.GERM[1],cols.GERM[3],T),
                 c("NEURO",cols.NEURO[1],cols.NEURO[3],F),
                 c("NDD",cols.NEURO[1],cols.NEURO[3],T),
                 c("PSYCH",cols.NEURO[1],cols.NEURO[3],F),
                 c("SOMA",cols.SOMA[1],cols.SOMA[3],F))
#Iterate and plot
lapply(plotVals,function(vals){
  df <- readData(vals[1],"E2","all")
  png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/",
            as.character(vals[1]),"_E2_all.mirrorManhattan.png",sep=""),
      height=1200,width=1500,res=300)
  mirrorManhattan(df,col.even=vals[2],col.odd=vals[3],ymax=20,yaxis=vals[4],adjusted=0.05/10420)
  dev.off()
})




