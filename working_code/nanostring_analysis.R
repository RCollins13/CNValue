#!/usr/bin/env Rscript

# Rare CNV Map
# Code copyright (c) 2018 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to perform analysis of Nanostring validation data for candidate triplosensitive genes


####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/SSC_LCL_nanostring_and_Zebrafish_dupValidation/rCNV_nanostring_analysis/"
PLOTDIR <- paste(WRKDIR,"/plots/",sep="")
if(!dir.exists(PLOTDIR)){
  dir.create(PLOTDIR)
}
sample.cols <- c("#C92836","#6A8CC6")
require(FactoMineR)
require(beeswarm)
require(vioplot)


#####################
#####Helper functions
#####################
#Scatterplot of raw vs. corrected expression values
correctionScatter <- function(gene){
  raw <- dat[,which(colnames(dat)==gene)]
  corrected <- corrected.expression.vals[,which(colnames(corrected.expression.vals)==gene)]
  plot(raw,corrected,lwd=2,
       xlab="",ylab="")
  mtext(1,text="Raw Expression (A.U.)",line=2.5)
  mtext(2,text="Corrected Expression (A.U.)",line=2.5)
  mtext(3,text=gene,font=2,line=0.5)
  abline(lm(corrected ~ raw),col="red")
  legend("topleft",bg=NA,col=NA,bty="n",
         legend=paste("R2 = ",round(cor(raw,corrected)^2,digits=3)))
}
#Single gene swarmplot, colored by status
singleGeneSwarm <- function(gene,expected.samples=NULL){
  #Get corrected expression values
  vals <- corrected.expression.vals[,which(colnames(corrected.expression.vals)==gene)]
  ylims <- range(vals)
  names(vals) <- dat$sample
  
  #Get expected sample indexes & colors
  sample.colors <- rep("gray85",times=length(vals))
  sample.outlines <- rep(NA,times=length(vals))
  if(is.null(expected.samples)){
    expected.samples <- which(unlist(lapply(expected.genes,function(genes){
      if(any(genes %in% gene)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })))
  }
  if(length(expected.samples)>0){
    expected.sample.colors <- dat$ASD[expected.samples]
    expected.sample.colors[which(expected.sample.colors=="Case")] <- sample.cols[1]
    expected.sample.colors[which(expected.sample.colors=="Control")] <- sample.cols[2]
    sample.colors[expected.samples] <- expected.sample.colors
    sample.outlines[expected.samples] <- "black"
  }
  
  #Prep plot area
  par(mar=c(1,4,2,0.5),bty="n")
  plot(x=c(0,1),y=c(ylims),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  
  #Add noise threshold
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=noise.thresh,
       col=adjustcolor("black",alpha=0.3),border=NA)
  abline(h=noise.thresh,lwd=2)
  
  #Add boxplot & dots
  if(length(expected.samples)>0){
    vioplot(vals[-expected.samples],add=T,wex=0.4,drawRect=F,col=NA,border="gray60",at=0.5)
  }else{
    vioplot(vals,add=T,wex=0.4,drawRect=F,col=NA,border="gray60",at=0.5)
  }
  beeswarm(vals,pch=21,add=T,at=0.5,corral="wrap",corralWidth=0.4,
           pwbg=sample.colors,lwd=2,pwcol=sample.outlines)
  
  #Dress up plot
  mtext(3,line=0.3,font=2,text=gene)
  legend("topright",bty="n",bg=NA,pch=19,col=c(sample.cols,"gray85"),
         legend=c("Affected CNV carrier",
                  "Unaffected non-carrier\nfamily control",
                  "Unrelated non-carriers"),
         cex=0.5)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),line=-0.4,tick=F,
       labels=axTicks(2),las=2,cex.axis=0.7)
  mtext(2,line=3,text="mRNA Expression (A.U.)")
}



######################
#####Read & clean data
######################
#Read master matrix
dat <- read.table(paste(WRKDIR,"rCNV_nanostring_count_matrix.wMetadata.txt",sep=""),header=T,sep="\t")
#Read gene list
genelist <- read.table(paste(WRKDIR,"nanostring_analysis_gene_list.txt",sep=""),header=T)
endogenous.genes <- genelist[which(genelist$class %in% c("Endogenous","Housekeeping")),1]
#Normalize all counts vs sum of positive/negative system controls and housekeeping genes
norm.genes.idx <- which(colnames(dat) %in% genelist$gene[which(genelist$class %in% c("Housekeeping","Negative","Positive"))])
norm.vect <- apply(dat[,norm.genes.idx],1,sum)
dat[,which(colnames(dat) %in% genelist$gene)] <- t(sapply(1:nrow(dat),function(i){
  i.vals <- as.numeric(dat[i,which(colnames(dat) %in% genelist$gene)])
  return(i.vals/norm.vect[i])
}))
#Compute IQR, mean, and sd of expression levels per sample
log.stdev.QC <- t(apply(dat[,which(colnames(dat) %in% endogenous.genes)],1,function(vals){
  logvals <- log10(as.numeric(vals))
  return(c(IQR(logvals),mean(logvals),sd(logvals)))
}))
#Exclude samples that are outliers
outlier.samples <- unique(unlist(sapply(1:3,function(i){
  vals <- log.stdev.QC[,i]
  outliers <- which(vals>quantile(vals,0.75)+(3*IQR(vals)) | 
                      vals<quantile(vals,0.25)-(3*IQR(vals)))
  return(outliers)
})))
dat <- dat[-outlier.samples,]
#Unlist expected genes per sample
expected.genes <- strsplit(dat$genes.all,split=",")
names(expected.genes) <- dat$sample


#############################################
#####Correct expression values for covariates
#############################################
#Compute PCA for expression data
pca <- PCA(dat[,which(colnames(dat) %in% genelist$gene)],graph=F,ncp=10)
#Prepare per-sample matrix for per-gene expression regression
sample.covariates <- data.frame()
for(i in 1:nrow(dat)){
  vals <- as.vector(dat[i,])
  names(vals) <- colnames(dat)
  
  #Dummy variables
  case.status <- length(which(vals$ASD %in% "Case"))
  sex <- length(which(vals$sex %in% "MALE"))
  european <- length(which(vals$ethnicity %in% "WHITE"))
  RNA.batch2 <- length(which(vals$batch.RNA %in% 2))
  RNA.batch3 <- length(which(vals$batch.RNA %in% 3))
  RNA.batch4 <- length(which(vals$batch.RNA %in% 4))
  cartridge.batch2 <- length(which(vals$batch.cartridge %in% 2))
  cartridge.batch3 <- length(which(vals$batch.cartridge %in% 3))
  cartridge.batch4 <- length(which(vals$batch.cartridge %in% 4))
  # QC.imaging <- length(which(vals$QC.imaging %in% 1))
  # QC.binding_density <- length(which(vals$QC.binding_density %in% 1))
  # QC.positive_control <- length(which(vals$QC.positive_control %in% 1))
  # QC.detection_limit <- length(which(vals$QC.detection_limit %in% 1))
  
  #Return vector
  sample.covariates <- rbind(sample.covariates,
                             c(vals$sample,case.status,sex,european,
                               RNA.batch2,RNA.batch3,RNA.batch4,
                               cartridge.batch2,cartridge.batch3,cartridge.batch4,
                               pca$ind$coord[i,]))
}
colnames(sample.covariates) <- c("sample","ASD","sex","european",
                                 "RNA.batch2","RNA.batch3","RNA.batch4",
                                 "cartridge.batch2","cartridge.batch3","cartridge.batch4",
                                 # "QC.imaging","QC.binding_density","QC.positive_control","QC.detection_limit",
                                 paste("PC",1:10,sep="."))
#Iterate per gene and calculate coefficients for each covariate
betas <- sapply(genelist$gene,function(gene){
  #Make data frame of expression values and covariates
  prefit.df <- cbind(dat[,which(colnames(dat)==gene)],
                     sample.covariates[,-1])
  colnames(prefit.df) <- c("expression",colnames(sample.covariates)[-1])
  prefit.df[,grep("PC.",colnames(prefit.df))] <- apply(prefit.df[,grep("PC.",colnames(prefit.df))],2,as.numeric)
  
  #Fit linear model
  fit <- lm(expression ~ ASD + sex + european 
            + RNA.batch2+ RNA.batch3 + RNA.batch4
            + cartridge.batch2 + cartridge.batch3 + cartridge.batch4
            + PC.1 + PC.2 + PC.3 + PC.4 + PC.5 + PC.6 + PC.7 + PC.8 + PC.9 + PC.10,
            data=prefit.df)
  
  #Return covariates
  return(fit$coefficients)
})
#Calculate corrected expression per sample per gene based on fit coefficients
corrected.expression.vals <- sapply(genelist$gene,function(gene){
  #Prep df
  prefit.df <- cbind(sample.covariates,
                     dat[,which(colnames(dat)==gene)])
  colnames(prefit.df) <- c(colnames(sample.covariates),"raw_expression")
  prefit.df[,-1] <- apply(prefit.df[,-1],2,as.numeric)
  
  #Apply fit per sample
  prefit.df$corrected_expression <- apply(prefit.df[,-1],1,function(vals){
    newval <- vals[length(vals)]+sum(betas[,which(colnames(betas)==gene)][-1]*vals[-length(vals)])
    newval <- max(c(0,newval))
    return(newval)
  })
})
#Generate per-gene plots of raw vs corrected expression
sapply(endogenous.genes,function(gene){
  pdf(paste(PLOTDIR,"/",gene,"_expression_correction_scatter.pdf",sep=""),
      height=4,width=4)
  par(mar=c(3.5,3.5,2,2))
  correctionScatter(gene)
  dev.off()
})


#####################################
#####Differential expression analysis
#####################################
#Calculate noise threshold based on negative controls - 5% FDR
noise.thresh <- quantile(as.vector(corrected.expression.vals[,grep("NEG_",colnames(corrected.expression.vals))]),0.95)
#Mean expression levels per gene
mean.gene.vals <- apply(corrected.expression.vals,2,mean)
#Get list of endogenous genes where mean is below noise threshold
length(which(mean.gene.vals[which(names(mean.gene.vals) %in% endogenous.genes)]<=noise.thresh))
#Get p-value per sample per endogenous gene
DE.p <- apply(corrected.expression.vals[,which(colnames(corrected.expression.vals) %in% endogenous.genes)],2,function(vals){
  zscores <- scale(vals,scale=T,center=T)
  pnorm(zscores[,1],lower.tail=F)
})
#FDR correct DE.p
DE.p.FDR <- matrix(p.adjust(DE.p,method="fdr"),byrow=F,nrow=nrow(dat))
colnames(DE.p.FDR) <- colnames(DE.p)
rownames(DE.p.FDR) <- dat$sample
#Get list of overexpressed genes per sample
DE.genes.per.sample.nom <- sapply(1:nrow(DE.p),function(i){
  paste(colnames(DE.p)[which(DE.p[i,]<=0.05)],collapse=",")
})
DE.genes.per.sample.FDR <- sapply(1:nrow(DE.p.FDR),function(i){
  paste(colnames(DE.p.FDR)[which(DE.p.FDR[i,]<=0.05)],collapse=",")
})
DE.genes.per.sample.bonf <- sapply(1:nrow(DE.p),function(i){
  paste(colnames(DE.p)[which(DE.p[i,]<=0.05/(length(endogenous.genes)*nrow(dat)))],collapse=",")
})
DE.genes.per.sample <- cbind("sample"=dat$sample,
                             "nominal"=DE.genes.per.sample.nom,
                             "FDR"=DE.genes.per.sample.FDR,
                             "Bonferroni"=DE.genes.per.sample.bonf)
#Write out results
write.table(DE.genes.per.sample,paste(WRKDIR,"Nanostring_DE_genes.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")


##########################
#####MASTER PLOTTING BLOCK
##########################
#Generate per-gene plots of expression per sample
sapply(endogenous.genes,function(gene){
  pdf(paste(PLOTDIR,"/",gene,"_gene_expression_distribution.pdf",sep=""),
      height=4,width=4)
  par(mar=c(3.5,3.5,2,2))
  singleGeneSwarm(gene)
  dev.off()
})
#Generate per-CNV plots of expression for all genes relevant to a given proband
sapply(which(dat$ASD=="Case"),function(i){
  genes <- unlist(strsplit(dat$genes.all[i],split=","))
  genes <- genes[which(genes %in% endogenous.genes)]
  genes.in.CNV <- unlist(strsplit(dat$genes.CNV[i],split=","))
  genes.in.CNV <- genes.in.CNV[which(genes.in.CNV %in% genes)]
  s.idx <- which(dat$family==dat$family[i])
  pdf(paste(PLOTDIR,"/",dat$family[i],"_CNV_interval_expression.pdf",sep=""),
      height=4,width=3+2*length(genes))
  par(mfrow=c(1,length(genes)),mar=c(3.5,3.5,2,2))
  sapply(genes,function(gene){
    singleGeneSwarm(gene,expected.samples=s.idx)
    if(gene %in% genes.in.CNV){
      mtext(1,line=0,text="IN CNV",font=2,col="red")
    }
  })
  dev.off()
})







