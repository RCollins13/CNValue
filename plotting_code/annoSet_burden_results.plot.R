#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmaps, dotplots, and volcano plots of results from annoSet burden testing

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Function to generate grid of dotplots for all tests
single.dotplot <- function(CNV,VF,filt,pheno,basedir){
  #Read data
  estimates <- read.table(paste(basedir,CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
  lower <- read.table(paste(basedir,CNV,"_",VF,"_",filt,".lowerCI.txt",sep=""),header=F)
  upper <- read.table(paste(basedir,CNV,"_",VF,"_",filt,".upperCI.txt",sep=""),header=F)
  pvals <- read.table(paste(basedir,CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)

  #Subset data to phenotype of interest and merge into single df
  results <- data.frame(estimates=estimates[,which(phenos %in% pheno)+1],
                        lower=lower[,which(phenos %in% pheno)+1],
                        upper=upper[,which(phenos %in% pheno)+1],
                        pvals=pvals[,which(phenos %in% pheno)+1])

  #Sort df by effect size point estimate & remove rows with NA estimates
  results <- results[order(results$estimates),]
  results <- results[which(!is.na(results$estimates)),]

  #Set color vector
  cols <- sapply(as.numeric(results$pvals),function(p){
    if(!is.na(p)){
      if(p<0.05/nrow(results)){
        return("red")
      }else{
        return("black")
      }
    }else{
      return("black")
    }
  })

  #Prepare plot & background
  par(mar=c(0.5,3,0.5,0.5))
  plot(x=c(0,nrow(results)+1),y=c(0,1.1*max(results$estimates,na.rm=T)),
       xaxs="i",xlab="",xaxt="n",yaxs="i",ylab="",yaxt="n",type="n")
  abline(h=axTicks(2),lwd=0.7,col="gray70")
  abline(h=1,lwd=2)

  #Plot point estimates and 95% CIs
  segments(x0=1:nrow(results),x1=1:nrow(results),
           y0=results$lower,y1=results$upper,
           lwd=0.5,col=cols)
  points(x=1:nrow(results),y=results$estimates,pch=18,cex=0.5,col=cols)

  #Add axes & labels
  mtext(1,line=-1,text=paste("Annotation Sets (n=",nrow(results),")",sep=""))
  axis(2,at=axTicks(2),las=2,line=-0.3,tick=F)
  axis(2,at=axTicks(2),las=2,labels=NA)
  mtext(2,line=1.6,text="CNV Burden")
  mtext(3,line=-1,font=2,text=pheno)

  #Add legend
  legend("topleft",legend=c("Effect Estimate","95% CI","Significant"),cex=0.5,
         pch=c(18,NA,18),lwd=c(NA,1,NA),col=c("black","black","red"),bg="white")
}

#####Code to generate all dotplot grids
par(mfrow=c(7,5))
sapply(phenos,function(pheno){
  single.dotplot("DEL","E2","all",pheno,
                 paste(WRKDIR,"plot_data/annoSet_burden_results/",sep=""))
})

