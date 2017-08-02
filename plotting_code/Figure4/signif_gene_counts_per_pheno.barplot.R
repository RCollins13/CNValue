#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate barplots of significant genes by pheno group (Fig 4)

###################
#####Set parameters
###################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("CTRL","GERM","NEURO","NDD","ASD","DD","ID","PSYCH","SCHIZ","BEHAV","SEIZ",
            "HYPO","UNK","NONN","DRU","GROW","SKIN","SKEL","HEAD","MUSC","HEART",
            "EE","EMI","CNCR","CLIV","CBRN","CBLD","CGST","CSKN","CRNL","CMSC",
            "CREP","CHNK","CBST","CLNG","CEND")
phenos.reorder <- c(1,c(1,12,
                    2,3,5,7,8,4,10,11,9,6,
                    13,18,15,20,17,14,19,21,16,22,
                    23,31,28,27,29,25,34,33,35,32,24,30,26)+1)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")
cols.allPhenos <- c(cols.CTRL[1],
                    cols.GERM[1],
                    rep(cols.NEURO[1],10),
                    cols.GERM[1],
                    rep(cols.SOMA[1],10),
                    rep(cols.CNCR[1],13))

##############
#####Read data
##############
df <- read.table(paste(WRKDIR,"plot_data/figure4/signif_genes_by_pheno.count.txt",sep=""),header=F)
names(df) <- c("pheno","DEL","BOTH","DUP")
df <- df[phenos.reorder,]
df$SUM <- df$DEL+df$BOTH+df$DUP
df.sumOrder <- order(df$SUM,decreasing=T)

########################
####Function for barplot
########################
plotBars <- function(df){
  #Determine max x-value
  xmax <- max(df$SUM)

  #Prepare plot area
  par(mar=c(0.5,2.5,1.5,0.5),
      bty="n")
  plot(x=c(-20,1.02*xmax),y=c(0,-nrow(df)),type="n",
       xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Plot gridlines
  abline(v=seq(0,500,25),col=cols.CTRL[4],lwd=0.5)
  abline(v=seq(0,500,50),col=cols.CTRL[3],lwd=0.75)
  abline(v=seq(0,500,100),col=cols.CTRL[2])

  #Plot rectangles
  rect(xleft=0,xright=df$DEL[df.sumOrder],
       ybottom=-(1:nrow(df))+0.25,
       ytop=-(1:nrow(df))+0.75,
       col="red",lwd=0.5)
  rect(xleft=df$DEL[df.sumOrder],
       xright=df$DEL[df.sumOrder]+df$BOTH[df.sumOrder],
       ybottom=-(1:nrow(df))+0.25,
       ytop=-(1:nrow(df))+0.75,
       col=cols.CTRL[2],lwd=0.5)
  rect(xleft=df$DEL[df.sumOrder]+df$BOTH[df.sumOrder],
       xright=df$SUM[df.sumOrder],
       ybottom=-(1:nrow(df))+0.25,
       ytop=-(1:nrow(df))+0.75,
       col="blue",lwd=0.5)
  rect(xleft=0,xright=df$SUM[df.sumOrder],
       ybottom=-(1:nrow(df))+0.25,
       ytop=-(1:nrow(df))+0.75,
       col=NA)

  #Plot y-axis rectanges
  rect(xleft=-20,xright=0,
       ybottom=-(1:nrow(df)),
       ytop=-(1:nrow(df))+1,
       col=cols.allPhenos[phenos.reorder][df.sumOrder])

  #Clean up plot lines
  rect(xleft=-20,xright=0,
       ybottom=-nrow(df),ytop=0,
       col=NA)

  #Add x-axis
  axis(3,at=seq(0,xmax,25),labels=NA,tck=-0.01,col=cols.CTRL[2])
  axis(3,at=seq(0,xmax,50),labels=NA,tck=-0.0175,col=cols.CTRL[1])
  axis(3,at=seq(0,xmax+100,100),labels=NA,tck=-0.0225)
  axis(3,at=seq(0,xmax+100,100),tick=F,
       labels=seq(0,xmax+100,100),line=-0.7)

  #Add y-axis
  sapply(df.sumOrder,function(i){
    axis(2,at=-i+0.5,tick=F,line=-0.8,las=2,cex.axis=0.75,
         labels=phenos[phenos.reorder][df.sumOrder][i])
  })
}

##############
#####Plot bars
##############
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/signif_gene_counts_per_pheno.barplots.pdf",sep=""),
    width=2.5,height=4.5)
plotBars(df)
dev.off()







