#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate density plots of significant element sizes (Fig 3)

###################
#####Set parameters
###################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)

######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

##############
#####Read data
##############
sizes <- lapply(list("GERM","NEURO","SOMA","CNCR"),function(pheno){
  dat <- read.table(paste(WRKDIR,"plot_data/NoncodingElementClassesFigure/element_size_distros/",
                          pheno,".element_sizes.txt",sep=""),header=F)[,1]
  dens <- sapply(seq(0,200000,1000),function(min){
    return(length(which(dat>min))/length(dat))
  })
})

#########
#####Plot
#########
#Prep plot area
par(mar=c(3,3,0.5,0.5))
plot(x=c(5,200),y=c(0,1),type="n",
     xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
#Add gridlines
abline(v=seq(0,200,12.5),col=cols.CTRL[4],lwd=0.5)
abline(v=seq(0,200,25),col=cols.CTRL[3],lwd=0.75)
abline(v=seq(0,200,50),col=cols.CTRL[2])
abline(h=seq(0,1,0.05),col=cols.CTRL[4],lwd=0.5)
abline(h=seq(0,1,0.1),col=cols.CTRL[3],lwd=0.75)
abline(h=seq(0,1,0.2),col=cols.CTRL[2])
#Plot lines
points(sizes[[1]],type="l",lwd=3,col=cols.GERM[1])
points(sizes[[2]],type="l",lwd=3,col=cols.NEURO[1])
points(sizes[[3]],type="l",lwd=3,col=cols.SOMA[1])
points(sizes[[4]],type="l",lwd=3,col=cols.CNCR[1])
#Add axes
axis(1,at=seq(0,200,25),labels=NA,tck=-0.01,col=cols.CTRL[1])
axis(1,at=seq(0,200,50),labels=NA,tck=-0.02)
axis(1,at=seq(0,200,50),labels=paste(seq(0,200,50),"kb",sep=""),tick=F,line=-0.5)
axis(2,at=seq(0,1,0.1),labels=NA,tck=-0.01,col=cols.CTRL[1])
axis(2,at=seq(0,1,0.2),labels=NA,tck=-0.02)
axis(2,at=seq(0,1,0.2),labels=paste(seq(0,100,20),"%",sep=""),tick=F,line=-0.5,las=2)



