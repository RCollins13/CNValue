#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate barplot of gene sets significantly enriched in CNV/DEL/DUP by phenotype

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

#####Load libraries
require(plotrix)

#####Load & clean data
#Read p-values
pvals <- lapply(c("CNV","DEL","DUP"),function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",CNV,
                          "_E2_all_exonic.pvals.txt",sep=""),header=F)
  dat[,-1] <- apply(dat[,-1],2,as.numeric)
  colnames(dat) <- c("geneSet",phenos)
  return(dat)
})
#Iterate over CNV types and count # signif phenotypes
signif <- lapply(pvals,function(dat){
  apply(dat[,-1],1,function(vals){
    return(length(which(vals<=0.05/263 & !is.na(vals))))
  })
})
#Bin # significant sets by 0, 1, 2-34, 35
signif.counts <- lapply(signif,function(vals){
  none <- length(which(vals==0))
  one <- length(which(vals==1))
  some <- length(which(vals>1 & vals<35))
  all <- length(which(vals==35))
  return(c(none,one,some,all))
})
#Count number of significant sets per pheno group
pheno.signif <- lapply(pvals,function(dat){
  apply(dat[,-1],2,function(vals){
    return(length(which(vals<=0.05/263 & !is.na(vals))))
  })
})
#Make master df of all counts for plotting on single axis
counts <- rbind(setNames(data.frame(signif.counts[[1]],
                                    signif.counts[[2]],
                                    signif.counts[[3]]),
                         c("CNV","DEL","DUP")),
                setNames(data.frame(pheno.signif[[1]],
                                    pheno.signif[[2]],
                                    pheno.signif[[3]]),
                         c("CNV","DEL","DUP")))
rownames(counts) <- c("None","One","Some","All",phenos)
#Adjust "none" category by -40 to fit on broken Y-axis (lim=100)
counts[1,] <- counts[1,]-40

#####Prepare plot area
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/geneSet_burden_count_bins.barplot.pdf",sep=""),
    width=2.7,height=4)
par(bty="n",mar=c(0.5,2.5,0.5,0.5))
plot(x=c(0,4),y=c(0,110),type="n",
     xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")

#####Draw gridlines
abline(h=seq(0,100,20),col=cols.CTRL[2])
abline(h=seq(10,90,20),col=cols.CTRL[3],lty=2)

#####Draw y-axis
axis(2,at=seq(0,90,10),col=cols.CTRL[1],
     labels=NA,tck=-0.015)
axis(2,at=seq(0,100,20),labels=NA,tick=-0.025)
axis(2,at=seq(0,100,20),tick=F,cex.axis=1.25,line=-0.3,
     labels=c(0,20,40,60,80,140),las=2)

#####Iterate over counts and plot bars
for(i in 1:4){
  rect(xleft=i-c(1,0.75,0.5)+0.125,
       xright=i-c(0.75,0.5,0.25)+0.125,
       ybottom=rep(0,3),ytop=as.numeric(counts[i,]),
       col=c("gray20","red","blue"))
}

#####Drop axis break
# rect(xleft=par("usr")[1],xright=par("usr")[2],
#      ybottom=93.5,ytop=94.5,
#      col="white",border=NA)
# abline(h=94,lty=2)
axis.break(2,breakpos=94,style="gap")

#####Drop border line
abline(h=0)

#####Close device
dev.off()

#####



