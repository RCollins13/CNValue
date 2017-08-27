#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate barplots of germline prevalance estimates for rCNV segments

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
#Mini helper function to read data
readData <- function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/figure2/prevalance/GERM_",
                          CNV,".prevalance.txt",sep=""),header=F)
  colnames(dat) <- c("pheno","hits","all")
  #Reorder
  dat <- dat[c(1,phenos.reorder+1)[1:23],]
  return(dat)
}
#Read data
DEL <- readData("DEL")
DUP <- readData("DUP")
#Merge data
CNV <- data.frame("pheno"=DEL[,1],
                  "DEL"=DEL[,2],
                  "DUP"=DUP[,2],
                  "CNV"=DEL[,2]+DUP[,2],
                  "N"=DEL[,3])
#Normalize data
CNV.norm <- CNV
CNV.norm[,2] <- CNV.norm[,2]/CNV.norm[,5]
CNV.norm[,3] <- CNV.norm[,3]/CNV.norm[,5]
CNV.norm[,4] <- CNV.norm[,4]/CNV.norm[,5]

###########################
#####Gather stats for paper
###########################
range(CNV.norm[-1,4])
binom.test(CNV[2,4],CNV[2,5],p=CNV.norm[1,4],alternative="greater")$p.value

###########################
#####Plot bars of CNV rates
###########################
#Prep plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure2/CNVsegments_incidence.barplots.pdf",sep=""),
    height=4.8,width=5)
par(mar=c(0.5,4,2.5,0.5),bty="n")
xmax <- max(ceiling(CNV.norm[,2:3]*200))/200
plot(x=c(-0.08*xmax,2.2*xmax),y=c(0,-nrow(CNV.norm)),
     type="n",xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")

#Add background elements
abline(v=c(seq(0,xmax,0.01),seq(1.2*xmax,2.2*xmax,0.01)),col=cols.CTRL[3])
rect(xleft=-0.08*xmax,xright=-0.01*xmax,
     ybottom=-1:-nrow(CNV),ytop=0:-(nrow(CNV)-1),
     col=c(cols.CTRL[1],cols.allPhenos[phenos.reorder]))
abline(v=c(0,1.2*xmax),lwd=2)

#Add axes
axis(3,at=seq(0,xmax,0.005),col=cols.CTRL[1],tck=-0.015,labels=NA)
axis(3,at=seq(0,xmax,0.01),tck=-0.025,labels=NA)
axis(3,at=seq(0,xmax,0.01),tick=F,line=-0.5,
     labels=paste(0:(100*xmax),"%",sep=""))
axis(3,at=seq(1.2*xmax,2.2*xmax,0.005),col=cols.CTRL[1],tck=-0.015,labels=NA)
axis(3,at=seq(1.2*xmax,2.2*xmax,0.01),tck=-0.025,labels=NA)
axis(3,at=seq(1.2*xmax,2.2*xmax,0.01),tick=F,line=-0.5,
     labels=paste(0:(100*xmax),"%",sep=""))
sapply(1:nrow(CNV),function(i){
  axis(2,at=-i+0.5,tick=F,line=-0.9,labels=c("CTRL",phenos[phenos.reorder])[i],las=2)
})

#Plot DEL rectangles
rect(xleft=0,xright=CNV.norm[,2],
     ybottom=(-1:-nrow(CNV))+0.2,
     ytop=(0:(-nrow(CNV)+1)-0.2),
     col="red")
rect(xleft=0,xright=CNV.norm[1,2],
     ybottom=((-1:-nrow(CNV))+0.2)[-1],
     ytop=((0:(-nrow(CNV)+1)-0.2))[-1],
     col=cols.CTRL[1])

#Plot DUP rectangles
rect(xleft=1.2*xmax,xright=(1.2*xmax)+CNV.norm[,3],
     ybottom=(-1:-nrow(CNV))+0.2,
     ytop=(0:(-nrow(CNV)+1)-0.2),
     col="blue")
rect(xleft=1.2*xmax,xright=(1.2*xmax)+CNV.norm[1,3],
     ybottom=((-1:-nrow(CNV))+0.2)[-1],
     ytop=((0:(-nrow(CNV)+1)-0.2))[-1],
     col=cols.CTRL[1])

# #Add reference lines
abline(v=CNV.norm[1,2])
axis(3,at=CNV.norm[1,2],lwd=2,labels=NA,col="red")
abline(v=(1.2*xmax)+CNV.norm[1,3])
axis(3,at=(1.2*xmax)+CNV.norm[1,3],lwd=2,labels=NA,col="blue")

#Close device
dev.off()




