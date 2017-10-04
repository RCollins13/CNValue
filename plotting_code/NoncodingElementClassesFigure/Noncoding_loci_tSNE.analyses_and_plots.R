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
points(tsne.110$Y[which(optics.res$cluster==31),],col="red",pch=19,cex=0.5)
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
pvals.all <- sapply(1:n.categories,calcEnrichments)
colnames(pvals.all) <- paste("c",1:ncol(pvals.all),sep="")
rownames(pvals.all) <- colnames(dat[,-c(1:4)])
write.table(pvals.all,"~/scratch/cluster_enrichments.txt",col.names=T,row.names=T,quote=F,sep="\t")

#####Generate barplots for feature enrichment p-values per cluster, and mini-tSNE showing location
#Helper function
plotPvalBars <- function(cidx,features,labels,
                         pmax=25,bonf=0.05/nrow(pvals.all)){
  #Reorder features and labels based on p-value
  pvals.plot <- sapply(1:length(features),function(i){
    return(-log10(pvals.all[which(rownames(pvals.all)==features[i]),cidx]))
  })
  features <- features[rev(order(pvals.plot))]
  labels <- labels[rev(order(pvals.plot))]
  #Prep plot area
  par(mar=c(0.5,15,2.5,0.5),bty="n")
  plot(x=c(0,1.02*pmax),y=c(0,length(features)),type="n",
       xaxs="i",xlab="",xaxt="n",yaxs="i",ylab="",yaxt="n")
  #Draw gridlines
  abline(v=seq(0,pmax,5),col=cols.CTRL[2])
  #Draw rectangles
  sapply(1:length(features),function(i){
    rect(xleft=0,xright=min(c(pmax,-log10(pvals.all[which(rownames(pvals.all)==features[i]),cidx]))),
         ybottom=length(features)-i+0.2,ytop=length(features)-i+0.8,
         col="red")
  })
  #Draw significance thresholds
  abline(v=-log10(c(0.05,bonf)),lty=2)
  #Add y-axis category labels
  axis(2,at=rev(1:length(features))-0.5,line=-0.8,tick=F,
       labels=labels,las=2)
  #Add top x-axis (p-vals)
  axis(3,at=c(seq(0,pmax,5)),labels=NA)
  axis(3,at=c(seq(0,pmax,5)),tick=F,line=-0.5,
       labels=c(seq(0,pmax-5,5),paste(">",pmax,sep="")))
  axis(3,at=-log10(c(0.05,bonf)),tick=F,line=-1,labels=c("N","B"),cex.axis=0.75,font=3)
  mtext(3,line=1.5,text="Enrichment P-value")
}
#Helper function to generate mini tSNE
minitsne <- function(cidx){
  par(mar=c(0.2,0.2,0.2,0.2),bty="n",bg="transparent")
  plot(tsne.110$Y,pch=19,cex=0.1,col=cols.CTRL[1],
       xaxt="n",yaxt="n",xlab="",ylab="")
  points(tsne.110$Y[which(optics.res$cluster==cidx),],col="red",pch=19,cex=0.2)
  axis(3,at=seq(par("usr")[1],par("usr")[2],
                by=(par("usr")[2]-par("usr")[1])/10),
       tck=0.02,labels=NA,lwd=1.5)
  axis(2,at=seq(par("usr")[3],par("usr")[4],
                by=(par("usr")[4]-par("usr")[3])/10),
       tck=0.02,labels=NA,lwd=1.5)
}
#Generate plots
#15: TBR cluster
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/mini_tSNE_c15.png",sep=""),
   width=400,height=400,res=300)
minitsne(15)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/feature_enrichment_bars_c15.pdf",sep=""),
    height=0.4+0.25*(7),width=6)
plotPvalBars(cidx=15,
             features=c("TBRs_All_Primary_Tissues","Strong_TFBS_CTCF","Strong_TFBS_SMC3",
                        "TFBS_ZNF143","TFBS_RAD21","Highly_conserved_ChromHMM_Weak_Transcribed_Elements"),
             labels=c("Tissue-Conserved TAD Boundaries","Strong CTCF Sites","Strong SMC3 Sites",
                      "ZNF143 Sites","RAD21 Sites","Weak Transcription"))
dev.off()
#5: Canonical enhancer cluster
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/mini_tSNE_c5.png",sep=""),
    width=400,height=400,res=300)
minitsne(5)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/feature_enrichment_bars_c5.pdf",sep=""),
    height=0.4+0.25*(7),width=6)
plotPvalBars(cidx=5,
             features=c("Evolutionarily_conserved_elements_GERP",
                        "Conserved_H3K27ac","Conserved_DHS","Conserved_FIREs","Conserved_Enhancers",
                        "Strong_TFBS_EP300","Strong_TFBS_CEBPB"),
             labels=c("Primary Sequence Conservation","Tissue-Conserved H3K27ac Peaks",
                      "Tissue-Conserved DHS","Tissue-Conserved FIREs","Tissue-Conserved Enhancers",
                      "Strong EP300 Sites","Strong CEBPB Sites"))
dev.off()
#1: Super enhancer cluster
png(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/mini_tSNE_c1.png",sep=""),
    width=400,height=400,res=300)
minitsne(1)
dev.off()
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/NoncodingElementClasses/feature_enrichment_bars_c1.pdf",sep=""),
    height=0.4+0.25*(7),width=6)
plotPvalBars(cidx=1,
             features=c("BRAIN_SuperEnhancer","HEART_SuperEnhancer",
                        "Strong_TFBS_All_TFs","CpG_Islands","Early_replicating_regions",
                        "Conserved_ChromHMM_Genic_Enhancers_Class_1","Highly_conserved_eQTLs"),
             labels=c("Brain Super Enhancers","Heart Super Enhancers","Strong TF Sites (70 TFs)","CpG Islands",
                      "Early Replication Timing","Tissue-Conserved Genic Enhancers","Tissue-Conserved eQTLs"))
dev.off()
#Inactive immune chromatin duplications
plotPvalBars(cidx=12,
             features=c("Late_replicating_regions","Ultraconserved_noncoding_elements_UCNEs",
                        "IMMUNE_CompartmentB","IMMUNE_ChromHMM_Heterochromatin",
                        "ENDOCRINE_ChromHMM_Quiescent","IMMUNE_DMRs","INT_DUP_significant"),
             labels=c("Late Replication Timing","Ultraconserved Noncoding Elements",
                      "Immune Inactive (B) Compartments","Immune Heterochromatin",
                      "Endocrine Quiescent Chromatin","Immune Differential Methylation",
                      "INT DUP Burden"))
#Bivalent brain enhancers in NDDs
plotPvalBars(cidx=2,
             features=c("NDD_DEL_significant","NDD_DUP_significant",
                        "BRAIN_ChromHMM_BivalentEnhancer","BRAIN_ChromHMM_BivalentPoisedTSS",
                        "BRAIN_ChromHMM_WeakRepressedPolycomb","BRAIN_eQTLs"),
             labels=c("NDD DEL Burden","NDD DUP Burden",
                      "Brain Bivalent Enhancer","Brain Poised TSS",
                      "Brain Weak Polycomb Repression","Brain eQTLs"))
#Heart/embryo-specific enhancers
plotPvalBars(cidx=31,
             features=c("EMBRYO_CompartmentA","EMBRYO_FIRE","HEART_ChromHMM_ActiveEnhancer1",
                        "HEART_H3K27ac","HEART_eQTLs"),
             labels=c("Embryonic Active (A) Compartments",
                      "Embryonic FIREs","Heart Active Enhancers","Heart H3K27ac Peaks","Heart eQTLs"))
#Brain-inactive, skin-active elements w/brain cancer dups
plotPvalBars(cidx=86,
             features=c("CBRN_DUP_significant","BRAIN_CompartmentB","SKIN_FIRE",
                        "SKIN_ChromHMM_TSSFlankUpstream","BRAIN_ChromHMM_WeakRepressedPolycomb",
                        "SKIN_Enhancer"),
             labels=c("CBRN DUP Burden","Brain Inactive (B) Compartments",
                      "Skin FIREs","Skin Upstream TSS Flanks","Brain Weak Polycomb Repression",
                      "Skin Enhancers"))







