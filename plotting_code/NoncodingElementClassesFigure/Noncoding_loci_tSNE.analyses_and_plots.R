#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Exploration of master significant noncoding locus matrix

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
tissues <- c("BRAIN","DIGESTIVE","EMBRYO","ENDOCRINE","FAT","HEADNECK",
             "HEART","IMMUNE","LUNG","MUSCLE","RENAL","REPRO","SKELETAL","SKIN")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(Rtsne)
require(dbscan)
require(MASS)

#####Read & clean master matrix
dat <- read.table(paste(WRKDIR,"plot_data/NoncodingElementClassesFigure/",
                        "master_noncoding_loci.annotated_matrix.bed.gz",sep=""),
                  header=T,comment.char="")
dat.scaled <- dat
dat.scaled[,-c(1:4)] <- apply(dat[,-c(1:4)],2,function(vals){
  if(max(vals,na.rm=T)>0){
    newVals <- vals/max(vals,na.rm=T)
  }else{
    newVals <- vals
  }
  return(newVals)
})
dat.binary <- dat
dat.binary[,-c(1:4)] <- apply(dat[,-c(1:4)],2,function(vals){
  newVals <- vals
  newVals[which(vals>0)] <- 1
  return(as.numeric(newVals))
})

#####Explore data
#Set coloring vectors
DEL.sig <- apply(dat[,grep("DEL",colnames(dat))],1,function(vals){
  if(any(vals>0)){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
DUP.sig <- apply(dat[,grep("DUP",colnames(dat))],1,function(vals){
  if(any(vals>0)){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
CNV.cols <- sapply(1:nrow(dat),function(i){
  if(DEL.sig[i]==T){
    if(DUP.sig[i]==T){
      return("darkorchid4")
    }else{
      return("red")
    }
  }else{
    if(DUP.sig[i]==T){
      return("blue")
    }else{
      return("black")
    }
  }
})
# #Heatmap - scaled
# png("~/scratch/heatmap_test.png",height=2000,width=2000)
# heatmap(dat.scaled[,-c(1:4)])
# dev.off()
# #Heatmap - binary
# png("~/scratch/heatmap_test.binary.png",height=2000,width=2000)
# heatmap(as.matrix(dat.binary[,-c(1:4)]))
# dev.off()
#PCA
# PCs <- prcomp(x=dat.binary[,-c(1:4)])
# plot(PCs$x[,1:2])
#tSNE
# tsne.50 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=50,verbose=TRUE,max_iter=1000)
# tsne.70 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=70,verbose=TRUE,max_iter=1000)
# tsne.90 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=90,verbose=TRUE,max_iter=1000)
# tsne.110 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=110,verbose=TRUE,max_iter=5000)
#Read tSNE data from saved session
load(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/noncoding_loci_tSNE.rdata",sep=""))
# tsne.130 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=130,verbose=TRUE,max_iter=1000)
# tsne.150 <- Rtsne(dat.scaled[,-c(1:4)],dims=2,perplexity=150,verbose=TRUE,max_iter=1000)

# tsne.110.bin <- Rtsne(dat.binary[,-c(1:4)],dims=2,perplexity=110,verbose=TRUE,max_iter=500)
# tsne.110.raw <- Rtsne(dat[,-c(1:4)],dims=2,perplexity=110,verbose=TRUE,max_iter=500)

# par(mfrow=c(2,3))
# lapply(list(tsne.50,tsne.70,tsne.90,tsne.110,tsne.130,tsne.150),function(ts){
#   plot(ts$Y,pch=19,cex=0.5)
#   points(ts$Y[which(apply(dat[,grep("Super",colnames(dat))],1,sum)>0),],col="red",pch=19,cex=0.5)
#   points(ts$Y[which(apply(dat[,grep("TBR",colnames(dat))],1,sum)>0),],col="blue",pch=19,cex=0.5)
# })

#OPTICS clustering on tSNE with perplexity=110
#Titrated clustering parameters until >90% of samples were assigned to clusters
# and clusters were an appropriate visual fit for tSNE
optics.res <- extractDBSCAN(optics(tsne.110$Y,eps=1.3,minPts=8),eps_cl=1.3)
optics.res
n.categories <- length(table(optics.res$cluster))-1
optics.cols.pal <- c("gray80",
                     sample(rainbow_hcl(n.categories,c=95,l=82,start=280,end=120),
                            n.categories,replace=F))

#####Master plot of tSNE with cluster colors
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/master_noncoding_tSNE.png",sep=""),
    width=1200,height=1200,res=300)
par(mar=c(0.5,0.5,0.5,0.5),bty="n")
plot(tsne.110$Y,pch=19,cex=0.3,
     xaxt="n",yaxt="n",xlab="",ylab="",
     col=optics.cols.pal[optics.res$cluster+1])
axis(3,at=seq(par("usr")[1],par("usr")[2],
              by=(par("usr")[2]-par("usr")[1])/10),
     tck=0.01,labels=NA,lwd=1.5)
# mtext(3,line=0.1,text="t-SNE 1")
axis(2,at=seq(par("usr")[3],par("usr")[4],
              by=(par("usr")[4]-par("usr")[3])/10),
     tck=0.01,labels=NA,lwd=1.5)
# mtext(2,line=0.1,text="t-SNE 2")
dev.off()

#####Mini tSNE colored by GERM/CNCR
GERM.hit <- which(apply(dat[,c(5:26,40:61)],1,sum)>0)
CNCR.hit <- which(apply(dat[,c(27:39,62:74)],1,sum)>0)
BOTH.hit <- GERM.hit[which(GERM.hit %in% CNCR.hit)]
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/mini_noncoding_tSNE.CNCR_vs_GERM.png",sep=""),
    width=800,height=800,res=300)
par(mar=c(0.5,0.5,0.5,0.5),bty="n")
plot(tsne.110$Y,pch=19,cex=0.25,type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")
points(tsne.110$Y[CNCR.hit,],col=cols.CNCR[1],pch=19,cex=0.3)
points(tsne.110$Y[GERM.hit,],col=cols.GERM[1],pch=19,cex=0.3)
points(tsne.110$Y[BOTH.hit,],col="#FF6A09",pch=19,cex=0.3)
axis(3,at=seq(par("usr")[1],par("usr")[2],
              by=(par("usr")[2]-par("usr")[1])/10),
     tck=0.02,labels=NA,lwd=1.5)
axis(2,at=seq(par("usr")[3],par("usr")[4],
              by=(par("usr")[4]-par("usr")[3])/10),
     tck=0.02,labels=NA,lwd=1.5)
dev.off()

#####Mini tSNE colored by DEL/DUP
DEL.hit <- which(apply(dat[,grep("DEL",colnames(dat))],1,sum)>0)
DUP.hit <- which(apply(dat[,grep("DUP",colnames(dat))],1,sum)>0)
BOTH.hit <- DEL.hit[which(DEL.hit %in% DUP.hit)]
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/mini_noncoding_tSNE.DEL_vs_DUP.png",sep=""),
    width=800,height=800,res=300)
par(mar=c(0.5,0.5,0.5,0.5),bty="n")
plot(tsne.110$Y,pch=19,cex=0.25,type="n",
     xaxt="n",yaxt="n",xlab="",ylab="")
points(tsne.110$Y[DEL.hit,],col="red",pch=19,cex=0.3)
points(tsne.110$Y[DUP.hit,],col="blue",pch=19,cex=0.3)
points(tsne.110$Y[BOTH.hit,],col=cols.CTRL[2],pch=19,cex=0.3)
axis(3,at=seq(par("usr")[1],par("usr")[2],
              by=(par("usr")[2]-par("usr")[1])/10),
     tck=0.02,labels=NA,lwd=1.5)
axis(2,at=seq(par("usr")[3],par("usr")[4],
              by=(par("usr")[4]-par("usr")[3])/10),
     tck=0.02,labels=NA,lwd=1.5)
dev.off()

sapply(1:length(table(optics.res$cluster)),function(i){
  kvals <- tsne.110$Y[which(optics.res$cluster==i),]
  # z <- kde2d(kvals[,1],kvals[,2],n=100)
  # contour(z,drawlabels=FALSE,nlevels=3,add=TRUE)

})
#text(tsne.110$Y,labels=optics.res$cluster,cex=0.5,col=optics.cols.pal[optics.res$cluster+1])

plot(tsne.110$Y,pch=19,cex=0.3,
     xaxt="n",yaxt="n",xlab="",ylab="")
points(tsne.110$Y[which(optics.res$cluster==18),],col="red",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("Super",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("TBR",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("FIRE",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("PMD",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("StrongH3K27ac",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("Ultraconserved",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("HAR",colnames(dat))],1,sum)>0),],col="black",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("DEL",colnames(dat))],1,sum)>0),],col="red",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("DUP",colnames(dat))],1,sum)>0),],col="blue",pch=19,cex=0.5)
# points(tsne.110$Y[which(apply(dat[,grep("DEL",colnames(dat))],1,sum)>0),],col="red",pch=19,cex=0.5)

#####Helper function to calculate binomial p-values across all annotations for a given cluster
calcEnrichments <- function(cluster){
  #Get indexes for cluster
  clustIdx <- which(optics.res$cluster==cluster)
  otherIdx <- which(optics.res$cluster!=cluster & optics.res$cluster>0)
  #Iterate over each feature
  pvals <- apply(dat.binary[,-c(1:4)],2,function(vals){
    return(binom.test(x=length(which(vals[clustIdx]>0)),
                      n=length(vals[clustIdx]),
                      p=length(which(vals[otherIdx]>0))/length(vals[otherIdx]),
                      alternative="greater")$p.value)
  })
  return(pvals)
}

#####Calculate binomial p-values for all clusters & write to file
pvals.all <- t(sapply(1:n.categories,calcEnrichments))
colnames(pvals.all) <- 1:ncol(pvals.all)
write.table(pvals.all,"~/scratch/cluster_enrichments.txt",col.names=T,row.names=T,quote=F,sep="\t")

#####Generate barplots for feature enrichment p-values per cluster
#Helper function
plotPvalBars <- function(cidx,)



