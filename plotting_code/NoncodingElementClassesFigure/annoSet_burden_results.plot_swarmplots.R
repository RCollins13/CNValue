#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate swarm plots for figure on noncoding annotation set-based rCNV burdens

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
require(beeswarm)

##########################################
#####Helper function to read OR data frame
##########################################
#Performs log2-transformation & overwrites all data as NA
readOR <- function(CNV,VF,filt){
  #Read data frame
  OR <- read.table(paste(WRKDIR,"plot_data/NoncodingElementClassesFigure/burdenMatrices/",
                         CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
  #log2-transform & clean
  OR[,-1] <- apply(OR[,-1],2,function(vals){
    vals <- log2(vals)
    vals[which(is.infinite(vals))] <- NA
    vals[which(is.nan(vals))] <- NA
    return(vals)
  })
  #Add column names
  colnames(OR) <- c("class",phenos)
  #Return
  return(OR)
}

###########################################################################
#####Helper function to subset data frame based on vector of class keywords
###########################################################################
keywordFilter.df <- function(df,keywords){
  indexes <- unique(unlist(sapply(keywords,function(keyword){
    grep(keyword,df$class,ignore.case=T)
  })))
  return(df[indexes,])
}

#############
#####Analysis
#############
#####Get mean odds ratios across all noncoding classes for GERM and CNCR
DEL.OR <- readOR("DEL","E2","haplosufficient")
DUP.OR <- readOR("DUP","E2","haplosufficient")
GERM.DEL.mean <- apply(DEL.OR[,2:23],1,mean,na.rm=T)
GERM.DUP.mean <- apply(DUP.OR[,2:23],1,mean,na.rm=T)
CNCR.DEL.mean <- apply(DEL.OR[,24:36],1,mean,na.rm=T)
CNCR.DUP.mean <- apply(DUP.OR[,24:36],1,mean,na.rm=T)
set.seed(1111)
sim.sd <- sd(c(GERM.DEL.mean,GERM.DUP.mean),na.rm=T)
CTRL.DEL.sim <- rnorm(mean=0,sd=sim.sd,n=length(GERM.DEL.mean))
#M-W U-test for DEL & DUP vs. simulations
2^mean(GERM.DEL.mean,na.rm=T)
format(wilcox.test(GERM.DEL.mean)$p.value,scientific=T)
2^mean(GERM.DUP.mean,na.rm=T)
format(wilcox.test(GERM.DUP.mean)$p.value,scientific=T)
2^mean(CNCR.DEL.mean,na.rm=T)
format(wilcox.test(CNCR.DEL.mean)$p.value,scientific=T)
2^mean(CNCR.DUP.mean,na.rm=T)
format(wilcox.test(CNCR.DUP.mean)$p.value,scientific=T)

#####Get mean ORs per phenotype for tissue-defined elements
#Iterate over all phenos
tissue.means <- sapply(2:36,function(i){
  #Calculate mean for all tissue-defined elements
  all.mean <- mean(c(keywordFilter.df(DEL.OR,tissues)[,i],
                     keywordFilter.df(DUP.OR,tissues)[,i]),na.rm=T)
  #Iterate over each tissue and calculate mean
  each.mean <- sapply(tissues,function(tissue){
    spec.means <- mean(c(keywordFilter.df(DEL.OR,tissue)[,i],
                         keywordFilter.df(DUP.OR,tissue)[,i]),na.rm=T)
  })
  #Normalize each vs all means
  norm.means <- (2^each.mean)/(2^all.mean)

})
colnames(tissue.means) <- phenos












# #####Load test data
# OR <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/DEL_E3_haplosufficient.effectSizes.txt",sep=""),header=F)
# colnames(OR) <- c("anno",phenos)
# OR[,-1] <- log2(OR[,-1])
# is.na(OR[,-1]) <- sapply(OR[,-1],is.infinite)

# #####Set dev params
# universe <- OR$NDD
# highlights <- list(c(grep("conserved",OR[,1]),grep("Conserved",OR[,1])),
#                    grep("Highly_conserved",OR[,1]),
#                    grep("BRAIN",OR[,1]),
#                    grep("FIRE",OR[,1]),
#                    grep("TBR",OR[,1]),
#                    grep("Strong",OR[,1]),
#                    grep("Super",OR[,1]))
# highlight.col <- cols.NEURO[1]
# ylim <- c(log2(1/3),log2(3))
# cex <- 0.3
# yaxis <- T

#####Function to plot a strip of swarms
swarmStrip <- function(universe,
                       highlights,
                       highlight.col,
                       cex,
                       ylim,
                       yaxis=T){
  #Prepare plot
  plot(x=c(-0.5,length(highlights)+0.5),y=ylim,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")

  #Y axis & label, if optioned
  if(yaxis==T){
    axis(2,at=log2(c(1/3,1/2,3/4,1,1.5,2,3)),labels=NA)
    axis(2,at=log2(c(1/3,1/2,3/4,1,1.5,2,3)),tick=F,
         line=-0.3,labels=c(0.33,0.5,0.75,1,1.5,2,3),las=2)
    mtext(2,text="Odds Ratio",line=2.5)
  }

  #Gridlines
  abline(h=log2(c(1/3,1/2,3/4,1,1.5,2,3)),col=cols.CTRL[3])
  abline(h=0,lwd=2)

  #Plot universe annotation set
  beeswarm(universe,add=T,at=0,method="center",corral="random",
           pch=19,col=highlight.col,cex=cex)
  sapply(1:length(highlights),function(i){
    beeswarm(universe,add=T,at=i,method="center",corral="random",
             pch=19,col=cols.CTRL[2],cex=cex)
  })

  #Add universe median & IQR
  abline(h=summary(universe)[3],lty=2,col="gray40")
  segments(x0=-0.15,x1=0.15,
           y0=summary(universe)[3],
           y1=summary(universe)[3],
           lwd=2)
  rect(xleft=-0.07,xright=0.07,
       ybottom=summary(universe)[2],
       ytop=summary(universe)[5])

  #Iterate over sets to be highlighted
  sapply(1:length(highlights),function(i){
    if(!all(is.na(universe[highlights[[i]]]))){
      #Highlight median & IQR
      segments(x0=i-0.15,x1=i+0.15,
               y0=summary(universe[highlights[[i]]])[3],
               y1=summary(universe[highlights[[i]]])[3],
               lwd=2)
      rect(xleft=i-0.07,xright=i+0.07,
           ybottom=summary(universe[highlights[[i]]])[2],
           ytop=summary(universe[highlights[[i]]])[5])
      #Highlight data

      beeswarm(universe[highlights[[i]]],add=T,at=i,method="center",corral="random",
               pch=19,col=highlight.col,cex=cex)
    }
  })
}

# #####Set dev parameters
# CNV <- "DUP"
# VF <- "E2"
# filt <- "haplosufficient"
# highlights <- list(c(grep("conserved",OR[,1]),grep("Conserved",OR[,1])),
#                    grep("Highly_conserved",OR[,1]),
#                    grep("BRAIN",OR[,1]),
#                    grep("FIRE",OR[,1]),
#                    grep("TBR",OR[,1]),
#                    grep("Strong",OR[,1]),
#                    grep("Super",OR[,1]))
# cex <- 0.25
# labels=c("Cons","HCons","Brain","FIRE","TBR","StrongFunc","SuperEnh")

#####Function to generate plots for average of all NEURO, SOMA, and CNCR
tripleStrip <- function(OR,highlights,cex,labels,
                        ylim=c(log2(1/3),log2(3))){
  #Clean data
  colnames(OR) <- c("anno",phenos)
  OR[,-1] <- log2(OR[,-1])
  is.na(OR[,-1]) <- sapply(OR[,-1],is.infinite)

  #Average across all NEURO, SOMA, and CNCR
  OR$NEURO.mean <- apply(OR[,2:13],1,mean,na.rm=T)
  OR$SOMA.mean <- apply(OR[,14:23],1,mean,na.rm=T)
  OR$CNCR.mean <- apply(OR[,24:36],1,mean,na.rm=T)

  #Prepare plot layout
  par(mfrow=c(3,1),mar=c(1.5,4,0,0),oma=c(0,0,0.5,0.5))

  #Plot swarm strips
  swarmStrip(OR$NEURO.mean,highlights,cols.NEURO[1],
             cex=cex,ylim=ylim)
  axis(1,at=0:length(highlights),tick=F,line=-0.7,labels=c("AllClasses",labels))
  swarmStrip(OR$SOMA.mean,highlights,cols.SOMA[1],
             cex=cex,ylim=ylim)
  axis(1,at=0:length(highlights),tick=F,line=-0.7,labels=c("AllClasses",labels))
  swarmStrip(OR$CNCR.mean,highlights,cols.CNCR[1],
             cex=cex,ylim=ylim)
  axis(1,at=0:length(highlights),tick=F,line=-0.7,labels=c("AllClasses",labels))
}

#####Plot noncoding triple strips
OR <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                       "DUP_E2_haplosufficient.effectSizes.txt",sep=""),header=F)
highlights <- list(c(128:133),                    #Sequence conservation
                   c(1:28),                       #Tissue conservation
                   c(29:56),                      #Strong tissue conservation
                   grep("BRAIN",OR[,1]),          #Brain-derived elements
                   grep("FIRE",OR[,1]),           #FIREs
                   grep("TBR",OR[,1]),            #TBRs
                   grep("Strong",OR[,1]),         #Strong biochemical elements
                   c(134:204),                    #TFBS
                   grep("Super",OR[,1]),          #Super enhancers
                   grep("eQTL",OR[,1]),           #eQTLs
                   grep("GWAS",OR[,1]),           #GWAS loci
                   grep("H3K27ac",OR[,1]),        #H3K27ac peaks
                   grep("DHS",OR[,1]),            #DHS peaks
                   grep("ActiveEnhancer",OR[,1])  #Active enhancers
)
labels <- c("SeqCons","TissCons","TissHCons","Brain","FIRE","TBR",
            "StrongBiochem","TFBS","SuperEnh","eQTLs","GWAS",
            "H3K27ac","DNAse","ActiveEnh")
#Plot
tripleStrip(OR,highlights,cex=0.22,labels)




