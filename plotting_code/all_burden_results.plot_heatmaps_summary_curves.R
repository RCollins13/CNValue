#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmaps and dotplots of results from annoSet and geneSet burden testing

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color palettes
require(RColorBrewer)
cols.germ <- colorRampPalette(c("white","#7B2AB3"))(41)
cols.neuro <- colorRampPalette(c("white","#00BFF4"))(41)
cols.soma <- colorRampPalette(c("white","#EC008D"))(41)
cols.cncr <- colorRampPalette(c("white","#FFCB00"))(41)








########################################################
#####Function to generate grid of dotplots for all tests
########################################################
single.dotplot <- function(estimates,lower,upper,pvals,
                           CNV,VF,filt,pheno){

  #Subset data to phenotype of interest and merge into single df
  results <- data.frame(estimates=estimates[,which(phenos %in% pheno)+1],
                        lower=lower[,which(phenos %in% pheno)+1],
                        upper=upper[,which(phenos %in% pheno)+1],
                        pvals=pvals[,which(phenos %in% pheno)+1])

  #Sort df by effect size point estimate & remove rows with Inf or NA estimates
  results <- results[order(results$estimates),]
  results[which(is.infinite(results$estimates)),1] <- NA
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
  plot(x=c(0,nrow(results)+1),y=c(0,1.1*quantile(results$estimates,probs=0.99,na.rm=T)),
       xaxs="i",xlab="",xaxt="n",yaxs="i",ylab="",yaxt="n",type="n")
  abline(h=axTicks(2),lwd=0.7,col="gray70")
  abline(h=1,lwd=2)

  #Plot point estimates and 95% CIs
  segments(x0=1:nrow(results),x1=1:nrow(results),
           y0=results$lower,y1=results$upper,
           lwd=0.5,col=cols)
  points(x=1:nrow(results),y=results$estimates,pch=18,cex=0.5,col=cols)

  #Add point for median effect
  points(x=nrow(results)/2,y=median(results$estimates),pch=23,bg="yellow",cex=1.5)
  text(x=nrow(results)/2,y=median(results$estimates),pos=3,
       labels=round(median(results$estimates),3),font=2)

  #Add axes & labels
  mtext(1,line=-1,text=paste("Annotation Sets (n=",nrow(results),")",sep=""))
  axis(2,at=axTicks(2),las=2,line=-0.3,tick=F)
  axis(2,at=axTicks(2),las=2,labels=NA)
  mtext(2,line=1.8,text="Odds Ratio",cex=0.85)
  mtext(3,line=-1.2,text=paste(pheno,": ",CNV," (",VF,", ",filt,")",sep=""))
}







#######################################
#####Code to generate all dotplot grids
#######################################
#annoSet burdens
#One pdf per phenotype
sapply(phenos,function(pheno){
  #Prep pdf
  pdf(paste(WRKDIR,"plots/annoSet_burden_results/",pheno,
            "_annoSet_burden_results.dotplots.pdf",sep=""),
      height=7,width=12)
  par(mfrow=c(4,3))
  #Pages: coding, haplosufficient, noncoding
  sapply(c("all","haplosufficient","noncoding"),function(filt){
    #Rows: E2, E3, E4, N1
    sapply(c("E2","E3","E4","N1"),function(VF){
      #Columns: CNV, DEL, DUP
      sapply(c("CNV","DEL","DUP"),function(CNV){
        #Read data
        estimates <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                      CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
        lower <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".lowerCI.txt",sep=""),header=F)
        upper <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".upperCI.txt",sep=""),header=F)
        pvals <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)
        #Make dotplot
        single.dotplot(estimates,lower,upper,pvals,CNV,VF,filt,pheno)
      })
    })
  })
  #Close pdf
  dev.off()
})
#geneSet burdens
#One pdf per phenotype
sapply(phenos,function(pheno){
  #Prep pdf
  pdf(paste(WRKDIR,"plots/geneSet_burden_results/",pheno,
            "_geneSet_burden_results.dotplots.pdf",sep=""),
      height=7,width=12)
  par(mfrow=c(4,3))
  #Pages: exonic, wholegene
  sapply(c("exonic","wholegene"),function(context){
    #Rows: E2, E3, E4, N1
    sapply(c("E2","E3","E4","N1"),function(VF){
      #Columns: CNV, DEL, DUP
      sapply(c("CNV","DEL","DUP"),function(CNV){
        #Read data
        estimates <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                      CNV,"_",VF,"_all_",context,".effectSizes.txt",sep=""),header=F)
        lower <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                  CNV,"_",VF,"_all_",context,".lowerCI.txt",sep=""),header=F)
        upper <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                  CNV,"_",VF,"_all_",context,".upperCI.txt",sep=""),header=F)
        pvals <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                  CNV,"_",VF,"_all_",context,".pvals.txt",sep=""),header=F)
        #Make dotplot
        single.dotplot(estimates,lower,upper,pvals,CNV,VF,context,pheno)
      })
    })
  })
  #Close pdf
  dev.off()
})
#allSet burdens
#One pdf per phenotype
sapply(phenos,function(pheno){
  #Prep pdf
  pdf(paste(WRKDIR,"plots/allSet_burden_results/",pheno,
            "_allSet_burden_results.dotplots.pdf",sep=""),
      height=7,width=12)
  par(mfrow=c(4,3))
  #Pages: coding, haplosufficient, noncoding
  sapply(c("all","haplosufficient","noncoding"),function(filt){
    #Rows: E2, E3, E4, N1
    sapply(c("E2","E3","E4","N1"),function(VF){
      #Columns: CNV, DEL, DUP
      sapply(c("CNV","DEL","DUP"),function(CNV){
        #Read data
        estimates <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                      CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
        lower <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".lowerCI.txt",sep=""),header=F)
        upper <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".upperCI.txt",sep=""),header=F)
        pvals <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                  CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)
        #Make dotplot
        single.dotplot(estimates,lower,upper,pvals,CNV,VF,filt,pheno)
      })
    })
  })
  #Close pdf
  dev.off()
})








############################################################################
#####Function to plot heatmaps of phenotypes (columns) vs annotations (rows)
############################################################################
single.heatmap <- function(estimates,pvals,pfilter=F,reorder=NULL){
  #Reorder data, if optioned
  if(!is.null(reorder)){
    neworder.terms <- read.table(reorder,header=F)[,1]
    new.estimates <- estimates[which(estimates[,1]==neworder.terms[1]),]
    for(i in 2:length(neworder.terms)){
      new.estimates <- rbind(new.estimates,estimates[which(estimates[,1]==neworder.terms[i]),],
                             deparse.level=0,make.row.names=F)
    }
    estimates <- new.estimates
    new.pvals <- pvals[which(pvals[,1]==neworder.terms[1]),]
    for(i in 2:length(neworder.terms)){
      new.pvals <- rbind(new.pvals,pvals[which(pvals[,1]==neworder.terms[i]),],
                         deparse.level=0,make.row.names=F)
    }
    pvals <- new.pvals
  }

  #Set all estimates to ~ [0,41]
  #0 denotes NA, should be greyed out
  #41 represents ≥ 4-fold enriched
  estimates[,-1] <- apply(estimates[,-1],2,function(col){
    col[which(col<1)] <- 0.1
    col[which(is.na(col))] <- 0
    col[which(is.infinite(col))] <- 4
    col[which(col>4)] <- 4
    return(round(10*col,0)+1)
  })

  #Prepare plot
  plot(x=c(0,ncol(estimates)-1),y=c(0,-nrow(estimates)),
       xaxs="i",xlab="",xaxt="n",yaxs="i",ylab="",yaxt="n",type="n")

  #Iterate over estimates and plot rectangles
  sapply(2:ncol(estimates),function(col){
    if(col < 4){
      rect(xleft=col-2,xright=col-1,
           ytop=0:-(nrow(estimates)-1),
           ybottom=-1:-nrow(estimates),
           border=NA,col=cols.germ[estimates[,col]])
    }else if(col >= 4 && col < 14){
      rect(xleft=col-2,xright=col-1,
           ytop=0:-(nrow(estimates)-1),
           ybottom=-1:-nrow(estimates),
           border=NA,col=cols.neuro[estimates[,col]])
    }else if (col >= 14 && col < 24){
      rect(xleft=col-2,xright=col-1,
           ytop=0:-(nrow(estimates)-1),
           ybottom=-1:-nrow(estimates),
           border=NA,col=cols.soma[estimates[,col]])
    }else{
      rect(xleft=col-2,xright=col-1,
           ytop=0:-(nrow(estimates)-1),
           ybottom=-1:-nrow(estimates),
           border=NA,col=cols.cncr[estimates[,col]])
    }
  })

  #Iterate over estimates and grey out non-significant results (if optioned)
  if(pfilter==T){
    sapply(2:ncol(estimates),function(col){
      blackout <- which(as.numeric(pvals[,col])>=0.05/nrow(pvals) | is.na(pvals[,col]))
      rect(xleft=col-2,xright=col-1,
           ytop=-(blackout-1),
           ybottom=-blackout,
           border=NA,col="white")
    })
  }

  #Add axes
  axis(3,at=0.5:(ncol(estimates)-1.5),line=-0.9,las=2,tick=F,labels=phenos,cex.axis=0.3)
  axis(2,at=-0.5:-(nrow(estimates)-0.5),line=-0.9,las=2,tick=F,labels=estimates[,1],cex.axis=0.3)
}








##############################
#####Code to plot all heatmaps
##############################
#annoSet heatmaps
#One heatmap per combination of filter, VF, and CNV
sapply(c("all","haplosufficient","noncoding"),function(filt){
  sapply(c("E2","E3","E4","N1"),function(VF){
    sapply(c("CNV","DEL","DUP"),function(CNV){
      #Read data
      estimates <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                    CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
      pvals <- read.table(paste(WRKDIR,"plot_data/annoSet_burden_results/",
                                CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)
      #Original sort order
      pdf(paste(WRKDIR,"plots/annoSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byTissue.AllResults.pdf",sep=""),
          height=35,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F)
      dev.off()
      pdf(paste(WRKDIR,"plots/annoSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byTissue.SignifResults.pdf",sep=""),
          height=35,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T)
      dev.off()
      #New sort order
      pdf(paste(WRKDIR,"plots/annoSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byElementClass.AllResults.pdf",sep=""),
          height=35,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_noncoding_annotations.alternative_sort.list",sep=""))
      dev.off()
      pdf(paste(WRKDIR,"plots/annoSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byElementClass.SignifResults.pdf",sep=""),
          height=35,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_noncoding_annotations.alternative_sort.list",sep=""))
      dev.off()
    })
  })
})
#geneSet heatmaps
#One pdf per context, VF, and CNV
sapply(c("exonic","wholegene"),function(context){
  sapply(c("E2","E3","E4","N1"),function(VF){
    sapply(c("CNV","DEL","DUP"),function(CNV){
      #Read data
      estimates <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                    CNV,"_",VF,"_all_",context,".effectSizes.txt",sep=""),header=F)
      pvals <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                                CNV,"_",VF,"_all_",context,".pvals.txt",sep=""),header=F)
      #Original sort order
      pdf(paste(WRKDIR,"plots/geneSet_burden_results/",
                CNV,"_",VF,"_all_",context,"_heatmap.byTissue.AllResults.pdf",sep=""),
          height=20,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F)
      dev.off()
      pdf(paste(WRKDIR,"plots/geneSet_burden_results/",
                CNV,"_",VF,"_all_",context,"_heatmap.byTissue.SignifResults.pdf",sep=""),
          height=20,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T)
      dev.off()
      #New sort order (mostly for expression levels)
      pdf(paste(WRKDIR,"plots/geneSet_burden_results/",
                CNV,"_",VF,"_all_",context,"_heatmap.byClass.AllResults.pdf",sep=""),
          height=20,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_gene_sets.alternative_sort.list",sep=""))
      dev.off()
      pdf(paste(WRKDIR,"plots/geneSet_burden_results/",
                CNV,"_",VF,"_all_",context,"_heatmap.byClass.SignifResults.pdf",sep=""),
          height=20,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_gene_sets.alternative_sort.list",sep=""))
      dev.off()
    })
  })
})
#allSet heatmaps
#One heatmap per combination of filter, VF, and CNV
sapply(c("all","haplosufficient","noncoding"),function(filt){
  sapply(c("E2","E3","E4","N1"),function(VF){
    sapply(c("CNV","DEL","DUP"),function(CNV){
      #Read data
      estimates <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                    CNV,"_",VF,"_",filt,".effectSizes.txt",sep=""),header=F)
      pvals <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                                CNV,"_",VF,"_",filt,".pvals.txt",sep=""),header=F)
      #Original sort order
      pdf(paste(WRKDIR,"plots/allSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byTissue.AllResults.pdf",sep=""),
          height=50,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F)
      dev.off()
      pdf(paste(WRKDIR,"plots/allSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byTissue.SignifResults.pdf",sep=""),
          height=50,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T)
      dev.off()
      #New sort order
      pdf(paste(WRKDIR,"plots/allSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byClass.AllResults.pdf",sep=""),
          height=50,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=F,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_all_annotations.alternative_sort.list",sep=""))
      dev.off()
      pdf(paste(WRKDIR,"plots/allSet_burden_results/",
                CNV,"_",VF,"_",filt,"_heatmap.byClass.SignifResults.pdf",sep=""),
          height=50,width=4)
      par(mar=c(0.5,8,2,0.5))
      single.heatmap(estimates,pvals,pfilter=T,
                     reorder=paste(WRKDIR,"rCNVmap/misc/master_all_annotations.alternative_sort.list",sep=""))
      dev.off()
    })
  })
})





