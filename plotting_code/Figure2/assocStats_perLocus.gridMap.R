#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate gridmap of all significant CNV loci and association statistics

####################################
#####Set parameters & load libraries
####################################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")
phenos.reorder <- c(1,12,
                        2,3,5,7,8,4,10,11,9,6,
                        13,18,15,20,17,14,19,21,16,22,
                        23,31,28,27,29,25,34,33,35,32,24,30,26)
require(plotrix)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")
cols.allPhenos <- c(cols.GERM[1],
                    rep(cols.NEURO[1],10),
                    cols.GERM[1],
                    rep(cols.SOMA[1],10),
                    rep(cols.CNCR[1],13))

##############
#####Read data
##############
#P-values
DEL.p <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DEL_pVal.bed",sep=""),header=F)
colnames(DEL.p) <- c("chr","start","end",phenos)
DEL.p <- DEL.p[c(1:3,phenos.reorder+3)]
DUP.p <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DUP_pVal.bed",sep=""),header=F)
colnames(DUP.p) <- c("chr","start","end",phenos)
DUP.p <- DUP.p[c(1:3,phenos.reorder+3)]
#Odds ratios
DEL.OR <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DEL_OR.bed",sep=""),header=F)
colnames(DEL.OR) <- c("chr","start","end",phenos)
DEL.OR <- DEL.OR[c(1:3,phenos.reorder+3)]
DUP.OR <- read.table(paste(WRKDIR,"plot_data/figure2/DEL_DUP_union.E2_all.signif.filtered.DUP_OR.bed",sep=""),header=F)
colnames(DUP.OR) <- c("chr","start","end",phenos)
DUP.OR <- DUP.OR[c(1:3,phenos.reorder+3)]

#Dev: take top 20 sites for test plot
DEL.p <- head(DEL.p,20)
DEL.OR <- head(DEL.OR,20)
DUP.p <- head(DUP.p,20)
DUP.OR <- head(DUP.OR,20)

#################################################
#####Helper function to plot vertical semicircles
#################################################
semicircles <- function(x,y,cols,radii){
  #Left semicircle
  if(radii[1]=="N"){
    if(cols[1]!=cols.CTRL[1]){
      text(x=x-0.75,y=y,labels="N",font=2,col=cols[1])
    }
  }else{
    floating.pie(xpos=x,ypos=y,x=c(1,1),start=pi/2,
                 col=c(cols[1],NA),border=NA,
                 radius=as.numeric(radii[1]))
  }
  #Right semicircle
  if(radii[1]=="N"){
    if(cols[2]!=cols.CTRL[1]){
      text(x=x+0.75,y=y,labels="N",font=2,col=cols[2])
    }
  }else{
    floating.pie(xpos=x,ypos=y,x=c(1,1),start=pi/2,
                 col=c(NA,cols[2]),border=NA,
                 radius=as.numeric(radii[2]))
  }
}

#######################################
#####Prepare plotting area & background
#######################################
#Set graphical parameters
scale=0.3
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/assoc_stats_per_locus.gridPlot.pdf",sep=""),
    width=scale*length(phenos),height=scale*nrow(DEL.p))
par(mar=c(0.5,9.5,4,0.5))
#Each cell is a 3x3 square
plot(x=c(0,3*length(phenos)),y=c(0,-3*nrow(DEL.p)),
     type="n",xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
#Add gridlines
abline(v=seq(0,3*length(phenos),3),h=seq(0,-3*nrow(DEL.p),-3),col=cols.CTRL[1])
#Add locus coordinates as left Y-axis labels
axis(2,at=seq(0,-3*nrow(DEL.p),-3),labels=NA,col=cols.CTRL[1])
sapply(1:nrow(DEL.p),function(i){
  axis(2,at=(-3*i)+1.5,tick=F,line=-0.8,las=2,cex.axis=0.75,
       labels=paste(DEL.p$chr[i],":",
                    prettyNum(DEL.p$start[i],big.mark=","),"-",
                    prettyNum(DEL.p$end[i],big.mark=","),sep=""))
})
#Add phenotype groups as X-axis labels
sapply(1:length(phenos),function(i){
  axis(3,at=(3*i)-1.5,tick=F,line=-0.8,las=2,labels=phenos[phenos.reorder[i]])
})

#####################################
#####Iterate over loci and plot stats
#####################################
sapply(1:nrow(DEL.p),function(i){
  #Iterate over phenotypes & plot deletion stats
  sapply(4:(length(phenos)+3),function(p){
    #Get colors
    piecols <- rep(cols.CTRL[1],2)
    if(DEL.p[i,p]<10^-8){
      piecols[1] <- "red"
    }
    if(DUP.p[i,p]<10^-8){
      piecols[2] <- "blue"
    }
    #Get odds ratios
    ORs <- c(DEL.OR[i,p],DUP.OR[i,p])
    ORs[which(is.nan(ORs))] <- 0
    ORs <- (log10(ORs)/2)+0.5
    ORs[which(ORs>1.5 & !is.infinite(ORs))] <- 1.5
    ORs[which(is.infinite(ORs))] <- "N"
    #Plot circle
    semicircles(x=(3*(p-3))-1.5,y=(-3*i)+1.5,
                cols=piecols,radii=ORs)
  })
})

##################
#####Clean up plot
##################
#Reapply gridlines
abline(v=seq(0,3*length(phenos),3),h=seq(0,-3*nrow(DEL.p),-3),col=cols.CTRL[1])
dev.off()








