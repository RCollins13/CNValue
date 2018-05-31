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
gene.cols <- c("Endogenous"="#3375B9",
               "Housekeeping"="#7ED8F3",
               "Negative"="#494A4D",
               "Positive"="#FF0000")
require(FactoMineR)
require(beeswarm)
require(vioplot)
require(psych)


#####################
#####Helper functions
#####################
#Scatterplot of raw vs. corrected expression values
correctionScatter <- function(gene){
  raw <- dat[,which(colnames(dat)==gene)]
  corrected <- corrected.expression.vals[,which(colnames(corrected.expression.vals)==gene)]
  plot(raw,corrected,lwd=2,
       xlab="",ylab="")
  mtext(1,text="Norm. Expression (A.U.)",line=2.5)
  mtext(2,text="Adj. Norm. Expression (A.U.)",line=2.5)
  mtext(3,text=gene,font=2,line=0.5)
  abline(lm(corrected ~ raw),col="red")
  legend("topleft",bg=NA,col=NA,bty="n",
         legend=paste("R2 = ",round(cor(raw,corrected)^2,digits=3)))
}
#Single gene swarmplot, colored by status
singleGeneSwarm <- function(gene,expected.samples=NULL,legend=T){
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
       col=adjustcolor("black",alpha=0.1),border=NA)
  abline(h=noise.thresh,lwd=2,lty=2)
  
  #Add boxplot & dots
  if(length(expected.samples)>0){
    vioplot(vals[-expected.samples],add=T,wex=0.4,drawRect=F,col=NA,border="gray60",at=0.5)
  }else{
    vioplot(vals,add=T,wex=0.4,drawRect=F,col=NA,border="gray60",at=0.5)
  }
  beeswarm(vals,pch=21,add=T,at=0.5,corral="wrap",corralWidth=0.4,
           pwbg=sample.colors,lwd=2,pwcol=sample.outlines)
  
  #Dress up plot
  if(gene %in% names(mean.gene.vals[which(mean.gene.vals<=noise.thresh)])){
    mtext(3,line=0.3,font=3,text=gene,col="gray60")
  }else{
    mtext(3,line=0.3,font=2,text=gene)
  }
  
  if(legend==T){
    legend("topright",bty="n",bg=NA,pch=19,col=c(sample.cols,"gray85"),
           legend=c("Affected CNV carrier",
                    "Unaffected non-carrier\nfamily control",
                    "Unrelated non-carriers"),
           cex=0.5)
  }
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),line=-0.4,tick=F,
       labels=axTicks(2),las=2,cex.axis=0.7)
  mtext(2,line=3,text="mRNA Expression (A.U.)")
}
#Plot single Gaddygram of expression values
gaddy <- function(vals,ymin=NULL,ymax=NULL,colors=NULL){
  #Prep plot area
  par(mar=c(5,3,3,1),bty="n")
  ylims <- round(as.numeric(quantile(as.matrix(vals),probs=c(0.025,0.995))))
  if(!is.null(ymin)){
    if(ymin<ylims[2]){
      ylims[1] <- ymin
    }
  }
  if(!is.null(ymax)){
    if(ymax>ylims[1]){
      ylims[2] <- ymax
    }
  }
  plot(x=c(0,ncol(vals)),y=ylims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  
  #Iterate over genes and plot each
  sapply(1:ncol(vals),function(i){
    #Get gene information
    gene <- colnames(vals)[i]
    gvals <- sort(as.numeric(vals[,i]))
    gvals[which(is.infinite(gvals))] <- NA
    if(is.null(colors)){
      col <- gene.cols[which(names(gene.cols)==genelist$class[which(genelist$gene==gene)])]
    }else{
      col <- colors[i]
    }
    
    #Plot points & line for median
    xpos <- seq(i-0.8,i-0.2,by=0.6/(nrow(vals)-1))
    points(x=xpos,y=gvals,pch=19,cex=0.3,col=col)
    segments(x0=i-0.8,x1=i-0.2,y0=mean(gvals,na.rm=T),y1=mean(gvals,na.rm=T),col=col,lwd=1.5)
    
    #Add gene label
    axis(1,at=i-0.5,tick=F,line=-0.8,las=2,labels=gene,font=3,cex.axis=0.7,col.axis=col)
  })
  
  #Add x-axis
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,labels=axTicks(2),las=2,cex.axis=0.7)
}
#Plot comparison of case/control number of DE genes
DE.gene.case.control.comp.plot <- function(case.count,control.count){
  #Format variables
  case.count <- as.numeric(case.count)
  control.count <- as.numeric(control.count)
  
  #Prep plot area
  par(mar=c(2,3.5,2,1))
  ylims <- c(0,max(c(case.count,control.count)+1))
  plot(x=c(0,2),y=ylims,type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")
  
  #Dress up plot
  axis(1,at=0.5,tick=F,line=-0.8,labels="Cases",col.axis=sample.cols[1],font=2)
  axis(1,at=1.5,tick=F,line=-0.8,labels="Controls",col.axis=sample.cols[2],font=2)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),labels=axTicks(2),las=2,line=-0.4,tick=F,cex.axis=0.8)
  mtext(2,line=2,text="DE Genes")
  
  #Plot data for cases
  if(any(case.count>0)){
    vioplot(case.count,add=T,wex=0.4,drawRect=F,col=NA,border=sample.cols[1],at=0.5)
  }
  boxplot(case.count,outline=F,lwd=1,staplewex=0,lty=1,boxwex=0.3,add=T,at=0.5,col=NA,xaxt="n",yaxt="n")
  beeswarm(case.count,pch=19,add=T,at=0.5,corral="wrap",corralWidth=0.4,
           col=sample.cols[1])
  
  #Plot data for controls
  if(any(control.count>0)){
    vioplot(control.count,add=T,wex=0.4,drawRect=F,col=NA,border=sample.cols[2],at=1.5)
  }
  boxplot(control.count,outline=F,lwd=1,staplewex=0,lty=1,boxwex=0.3,add=T,at=1.5,col=NA,xaxt="n",yaxt="n")
  beeswarm(control.count,pch=19,add=T,at=1.5,corral="wrap",corralWidth=0.4,
           col=sample.cols[2])
  
  #Add p-value
  text(x=1,y=par("usr")[4],pos=1,
       labels=paste("P = ",round(suppressWarnings(wilcox.test(case.count,control.count))$p.value,digits=4),sep=""))
}


######################
#####Read & clean data
######################
#Read master matrix
dat <- read.table(paste(WRKDIR,"rCNV_nanostring_count_matrix.wMetadata.txt",sep=""),header=T,sep="\t")
#Read gene list
genelist <- read.table(paste(WRKDIR,"nanostring_analysis_gene_list.txt",sep=""),header=T)
endogenous.genes <- genelist[which(genelist$class %in% c("Endogenous","Housekeeping")),1]
#Normalize all counts vs geometric mean of just housekeeping genes
norm.genes.idx <- which(colnames(dat) %in% genelist$gene[which(genelist$class %in% c("Housekeeping"))])
norm.vect <- apply(dat[,norm.genes.idx],1,geometric.mean,na.rm=T)
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
#Read list of genes ordered by chr & pos
genome.ordered.genes.all <- read.table(paste(WRKDIR,"gene_symbols_ordered_by_chr_pos.hg19.txt",sep=""),
                                       header=F)[,1]


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


##################################
#####Expression distribution plots
##################################
#####Gaddygram of all expression values per gene
#Calculate mean expression level per gene
mean.gene.vals <- apply(corrected.expression.vals,2,mean)
mean.gene.vals.log <- log10(mean.gene.vals)
mean.gene.vals.log.sort <- sort(mean.gene.vals.log)
#Sort gene expression matrix based on mean expression value
corrected.expression.vals.sort <- corrected.expression.vals[,order(mean.gene.vals)]
#Get number of genes per magnitude range
gaddy.range <- floor(min(mean.gene.vals.log)):ceiling(max(mean.gene.vals.log))
gaddy.range.table <- sapply(gaddy.range,function(i){
  length(which(mean.gene.vals.log>i & mean.gene.vals.log<=(i+1)))
})
names(gaddy.range.table) <- gaddy.range
gaddy.range <- gaddy.range.table[which(gaddy.range.table>0)]
if(gaddy.range[length(gaddy.range)]==1){
  gaddy.range <- gaddy.range[-length(gaddy.range)]
  gaddy.range[length(gaddy.range)] <- gaddy.range[length(gaddy.range)]+1
}
gaddy.range.table <- data.frame(c(1,cumsum(gaddy.range[-length(gaddy.range)])+1),
                                cumsum(gaddy.range))
#Plot gaddygram panels
pdf(paste(PLOTDIR,"/Gaddygram.all_genes.multipanel.pdf",sep=""),
    height=3,width=12)
layout(matrix(1:length(gaddy.range),nrow=1,byrow=T),
       widths=8+as.numeric(gaddy.range))
sapply(1:length(gaddy.range),function(i){
  gaddy(vals=corrected.expression.vals.sort[,gaddy.range.table[i,1]:gaddy.range.table[i,2]]/10^(as.numeric(names(gaddy.range)[i])-1))
  if(i==1){
    mtext(2,line=1.8,text="Normalized Expression (A.U.)",cex=0.8)
  }
  axis(3,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA,line=0.3)
  mtext(3,line=0.5,text=substitute("x10" ^X, list(X=as.numeric(names(gaddy.range)[i])-1)),cex=0.8)
})
dev.off()

#####Gaddygram of housekeeping genes
pdf(paste(PLOTDIR,"/Gaddygram.housekeeping.pdf",sep=""),
    height=3,width=5)
gaddy(vals=corrected.expression.vals.sort[,which(colnames(corrected.expression.vals.sort) %in% genelist$gene[which(genelist$class=="Housekeeping")])],
      ymin=0,ymax=0.6)
mtext(2,line=1.8,text="Normalized Expression (A.U.)",cex=0.8)
# axis(3,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA,line=0.3)
# mtext(3,line=0.5,text=substitute("x10" ^X, list(X=-1)),cex=0.8)
dev.off()

#####Gaddygram of positive control genes
pdf(paste(PLOTDIR,"/Gaddygram.positive_control.pdf",sep=""),
    height=4,width=6)
gaddy(vals=log2(corrected.expression.vals.sort[,which(colnames(corrected.expression.vals.sort) %in% genelist$gene[which(genelist$class=="Positive")])]))
mtext(2,line=1.8,text="log2 Normalized Expression (A.U.)",cex=0.8)
dev.off()

#####Gaddygram of negative control genes & noise threshold
#Calculate noise threshold based on negative controls - 5% FDR
noise.thresh <- quantile(as.vector(corrected.expression.vals[,grep("NEG_",colnames(corrected.expression.vals))]),0.95)
#Get list of endogenous genes where mean is below or at noise threshold
length(which(mean.gene.vals[which(names(mean.gene.vals) %in% endogenous.genes)]<=noise.thresh))
pdf(paste(PLOTDIR,"/Gaddygram.negative_control.pdf",sep=""),
    height=4,width=6)
gaddy(vals=log2(corrected.expression.vals.sort[,which(colnames(corrected.expression.vals.sort) %in% genelist$gene[which(genelist$class=="Negative")])]))
rect(xleft=par("usr")[1],xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=log2(noise.thresh),
     border=NA,col=adjustcolor("black",alpha=0.2))
abline(h=log2(noise.thresh),lty=2)
mtext(2,line=1.8,text="log2 Normalized Expression (A.U.)",cex=0.8)
dev.off()

#####Gaddygram of lowest 43 endogenous genes w/noise threshold
gaddy.noise.dat <- corrected.expression.vals.sort[,which(colnames(corrected.expression.vals.sort) %in% endogenous.genes)][,1:40]
gaddy.noise.cols <- rep("#4DAC26",times=ncol(gaddy.noise.dat))
gaddy.noise.cols[which(apply(gaddy.noise.dat,2,mean)<=noise.thresh)] <- "#D01C8B"
pdf(paste(PLOTDIR,"/Gaddygram.noise_threshold_applied.pdf",sep=""),
    height=3.5,width=7)
gaddy(vals=(10^3)*gaddy.noise.dat,colors=gaddy.noise.cols)
rect(xleft=par("usr")[1],xright=par("usr")[2],
     ybottom=par("usr")[3],ytop=(10^3)*noise.thresh,
     border=NA,col=adjustcolor("black",alpha=0.2))
abline(h=(10^3)*noise.thresh,lty=2)
axis(3,at=c(par("usr")[1],par("usr")[2]),tck=0,labels=NA,line=0.3)
mtext(3,line=0.5,text=substitute("x10" ^X, list(X=-3)),cex=0.8)
mtext(2,line=1.8,text="Normalized Expression (A.U.)",cex=0.8)
dev.off()


#####################################
#####Differential expression analysis
#####################################
#Get p-value per sample per endogenous gene
DE.p <- apply(corrected.expression.vals[,which(colnames(corrected.expression.vals) %in% endogenous.genes)],2,function(vals){
  zscores <- scale(vals,scale=T,center=T)
  pnorm(zscores[,1],lower.tail=F)
})
#Exclude genes that aren't expressed above noise threshold
DE.p <- DE.p[,which(!(colnames(DE.p) %in% names(which(names(mean.gene.vals) %in% endogenous.genes & mean.gene.vals<=noise.thresh))))]
#FDR correct DE.p
DE.p.FDR <- matrix(p.adjust(DE.p,method="fdr"),byrow=F,nrow=nrow(dat))
colnames(DE.p.FDR) <- colnames(DE.p)
rownames(DE.p.FDR) <- dat$sample
#Get list of overexpressed genes per sample
DE.genes.per.sample.nom <- sapply(1:nrow(DE.p),function(i){
  paste(colnames(DE.p)[which(DE.p[i,]<=0.05)],collapse=",")
})
DE.genes.per.sample.nom.ovr <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.nom[i],split=",")),
                  unlist(strsplit(dat$genes.all[i],split=","))),collapse=",")
})
DE.genes.per.sample.nom.ovr.key <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.nom[i],split=",")),
                  unlist(strsplit(dat$genes.key[i],split=","))),collapse=",")
})
DE.genes.per.sample.FDR <- sapply(1:nrow(DE.p.FDR),function(i){
  paste(colnames(DE.p.FDR)[which(DE.p.FDR[i,]<=0.05)],collapse=",")
})
DE.genes.per.sample.FDR.ovr <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.FDR[i],split=",")),
                  unlist(strsplit(dat$genes.all[i],split=","))),collapse=",")
})
DE.genes.per.sample.FDR.ovr.key <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.FDR[i],split=",")),
                  unlist(strsplit(dat$genes.key[i],split=","))),collapse=",")
})
DE.genes.per.sample.bonf <- sapply(1:nrow(DE.p),function(i){
  paste(colnames(DE.p)[which(DE.p[i,]<=0.05/(ncol(DE.p)*nrow(dat)))],collapse=",")
})
DE.genes.per.sample.bonf.ovr <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.bonf[i],split=",")),
                  unlist(strsplit(dat$genes.all[i],split=","))),collapse=",")
})
DE.genes.per.sample.bonf.ovr.key <- sapply(1:nrow(DE.p),function(i){
  paste(intersect(unlist(strsplit(DE.genes.per.sample.bonf[i],split=",")),
                  unlist(strsplit(dat$genes.key[i],split=","))),collapse=",")
})
DE.genes.per.sample <- cbind("sample"=dat$sample,
                             "family"=dat$family,
                             "ASD"=dat$ASD,
                             "all.genes.expected"=dat$genes.all,
                             "all.genes.key"=dat$genes.key,
                             "nominal"=DE.genes.per.sample.nom,
                             "nominal.expected"=DE.genes.per.sample.nom.ovr,
                             "nominal.key"=DE.genes.per.sample.nom.ovr.key,
                             "FDR"=DE.genes.per.sample.FDR,
                             "FDR.expected"=DE.genes.per.sample.FDR.ovr,
                             "FDR.key"=DE.genes.per.sample.FDR.ovr.key,
                             "Bonferroni"=DE.genes.per.sample.bonf,
                             "Bonferroni.expected"=DE.genes.per.sample.bonf.ovr,
                             "Bonferroni.key"=DE.genes.per.sample.bonf.ovr.key)
#Iterate over columns and restrict to genes with mean expression above the noise threshold
DE.genes.per.sample[,-c(1:3)] <- apply(DE.genes.per.sample[,-c(1:3)],2,function(vals){
  sapply(vals,function(genes){
    genes.v <- unlist(strsplit(genes,split=","))
    genes.filt <- names(which(names(mean.gene.vals) %in% genes.v & mean.gene.vals>noise.thresh))
    return(paste(sort(genes.filt),collapse=","))
  })
})
#Write out results
write.table(DE.genes.per.sample,paste(WRKDIR,"Nanostring_DE_genes.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

#####DE gene case vs control plot (six panels)
#Get count of nom/fdr/bonf sig DE genes per sample
DE.gene.counts.per.sample <- apply(DE.genes.per.sample[,-c(1:5)],2,function(strings){
  as.numeric(sapply(strings,function(str){
    length(unlist(strsplit(str,split=",")))
  }))
})
DE.gene.counts.per.sample <- as.data.frame(cbind(DE.genes.per.sample[,1:5],
                                                 apply(DE.gene.counts.per.sample,2,as.numeric)))
#Prep plot area
pdf(paste(PLOTDIR,"/DE_gene_counts.case_control.six_panel.pdf",sep=""),
    height=5,width=9)
par(mfrow=c(2,3))
DE.gene.case.control.comp.plot(case.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$nominal,
                               control.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$nominal)
mtext(3,line=0.5,text="Nominal",font=2)
DE.gene.case.control.comp.plot(case.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$FDR,
                               control.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$FDR)
mtext(3,line=0.5,text="FDR",font=2)
DE.gene.case.control.comp.plot(case.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$Bonferroni,
                               control.count=DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$Bonferroni)
mtext(3,line=0.5,text="Bonferroni",font=2)
DE.gene.case.control.comp.plot(case.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$nominal)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$nominal.expected),
                               control.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$nominal)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$nominal.expected))
DE.gene.case.control.comp.plot(case.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$FDR)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$FDR.expected),
                               control.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$FDR)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$FDR.expected))
DE.gene.case.control.comp.plot(case.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$Bonferroni)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Case"),]$Bonferroni.expected),
                               control.count=as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$Bonferroni)-
                                 as.numeric(DE.gene.counts.per.sample[which(DE.gene.counts.per.sample$ASD=="Control"),]$Bonferroni.expected))
dev.off()


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
  genes <- c(genome.ordered.genes.all[which(genome.ordered.genes.all %in% genes)],
             genes[which(!(genes %in% genome.ordered.genes.all))])
  genes.in.CNV <- unlist(strsplit(dat$genes.CNV[i],split=","))
  genes.in.CNV <- genes.in.CNV[which(genes.in.CNV %in% genes)]
  s.idx <- which(dat$family==dat$family[i])
  pdf(paste(PLOTDIR,"/",dat$family[i],"_CNV_interval_expression.pdf",sep=""),
      height=3,width=2+1.5*length(genes))
  par(mfrow=c(1,length(genes)),mar=c(3.5,3.5,2,2))
  sapply(genes,function(gene){
    singleGeneSwarm(gene,expected.samples=s.idx,legend=F)
    if(gene %in% genes.in.CNV){
      mtext(1,line=0,text="IN CNV",font=2,col="red")
    }
  })
  dev.off()
})







