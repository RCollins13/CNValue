#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#annoBurden: permutation testing for enrichment of annotations
# versus a specific subset of bins (e.g. those passing TBRden permutation
# for a significant case:control enrichment of CNVs)

annoBurden <- function(anno,             #Path to TBRden_pileup.sh output for the control group. BED format required
                       testBins,         #Path to TBRden_pileup.sh output for the comparison group. BED format required
                       perm=10000,       #Number of permutations
                       measure="mean",   #Measurement to use (options: mean, count; count computes only 0 | non-0)
                       alt="greater",    #Alternative hypothesis: test bins are (greater|less) than random. Two-tailed test: two.sided
                       column=5,         #Specify which column in df anno contains the annotation variable
                       plot=F,           #Option to generate barplot on output
                       OUTDIR=NULL,      #Output directory for writing results
                       prefix="TBRden"   #Prefix to be appended to results and plot (if optioned)
){
  #Sanity check input files
  if(!(file.exists(anno))){
    stop(paste("Input file ",anno," does not exist.",sep=""))
  }
  if(!(file.exists(testBins))){
    stop(paste("Input file ",testBins," does not exist.",sep=""))
  }

  #Sanity check arguments
  if(!(measure %in% c("mean","count"))){
    stop("Argument 'measure' must be either 'mean' or 'count'")
  }

  #Load library
  require(plotrix)

  #Read data
  anno <- read.table(anno,header=F,sep="\t")
  bins <- read.table(testBins,header=F,sep="\t")[,1:3]

  #Get annotation values corresponding to test bins
  anno.t <- merge(anno,bins,by=1:3,sort=F)

  #Sanity check to ensure the number of test bins matches number of test annotation bins
  if(nrow(anno.t) != nrow(bins)){
    warning(paste("Annotation bins do not match all test bins; only ",nrow(anno.t),
                  " of ",nrow(bins)," matched",sep=""))
  }

  #Slice annotation vector for test bins
  anno.t <- anno.t[,column]

  #Compute baseline statistics
  if(measure=="mean"){
    #Run t-test against background of all possible bins
    base <- t.test(anno.t,anno[,column],alternative=alt)
  }else{
    #Run binomial test against background of all possible bins
    x <- data.frame("nonzero"=c(length(which(anno.t>0)),
                                length(which(anno[,column]>0))),
                    "zero"=c(length(which(anno.t==0)),
                             length(which(anno[,column]==0))))
    base <- binom.test(x=x[1,1],
                       n=sum(x[1,]),
                       p=x[2,1]/sum(x[2,]),
                       alternative=alt)
  }

  #Apply over number of permutations
  res <- as.data.frame(t(sapply(1:perm,function(i){
    #Sample set of bins of same length as anno.t
    sim <- anno[sample(x=1:nrow(anno),size=length(anno.t),replace=F),column]

    #Perform either two-sample t-test or binomial test depending on specifications
    if(measure=="mean"){
      #Run t-test
      res <- t.test(sim,anno[,column],alternative=alt)

      #Return p-value & estimated mean
      return(c(res$p.value,res$estimate[1]))

    }else{
      #Build df of zero/nonzero bins
      x <- data.frame("nonzero"=c(length(which(sim>0)),
                                  length(which(anno[,column]>0))),
                      "zero"=c(length(which(sim==0)),
                               length(which(anno[,column]==0))))

      #Run binomial test
      res <- binom.test(x=x[1,1],
                        n=sum(x[1,]),
                        p=x[2,1]/sum(x[2,]),
                        alternative=alt)

      #Return p-value & estimated Bernoulli success rate
      return(c(res$p.value,res$estimate))
    }
  })))

  #Calculate empirical p-value & sampling distribution parameters
  p.emp <- (length(which(res[,1]<=base$p.value))+1)/(perm+1)
  est <- mean(res[,2])
  std.err <- std.error(res[,2],na.rm=T)

  #Barplot (if optioned)
  if(plot==T){
    if(measure=="mean"){
      #Prep barplot
      pdf(paste(OUTDIR,"/",prefix,".permutations.pdf",sep=""),
          height=6,width=2)
      par(mar=c(0.6,4.1,0.6,0.6))
      plot(x=c(0,2),y=c(0,1.2*max(est,base$estimate)),type="n",
           xaxs="i",yaxs="i",xlab="",xaxt="n",
           ylab="Elements per Bin (Mean)")
      rect(xleft=par("usr")[1],xright=par("usr")[2],
           ybottom=par("usr")[3],ytop=par("usr")[4],
           border=NA,col="gray95")
      yticks=seq(0,ceiling(par("usr")[4]),by=ceiling(par("usr")[4]/7))
      abline(h=yticks,lty=3)
      #Plot bars
      rect(xleft=0.2,xright=1,ybottom=0,ytop=base$estimate[1],
           col="darkorange1",border="darkorange3")
      rect(xleft=1,xright=1.8,ybottom=0,ytop=est,
           col="dodgerblue3",border="dodgerblue4")
      segments(x0=c(1.2,1.2,1.4),x1=c(1.6,1.6,1.4),
               y0=c(est-2*std.err,est-2*std.err,est+2*std.err),
               y1=c(est-2*std.err,est+2*std.err,est+2*std.err))
      #Add labels
      text(x=c(0.6,1.4),y=c(0,0),pos=3,labels=c("O","P"),font=2)
      dev.off()

    }else{
      #Prep strip plot
      pdf(paste(OUTDIR,"/",prefix,".permutations.pdf",sep=""),
          height=4,width=15)
      h <- hist(res[,2],breaks=0:100/100,plot=F)
      par(mar=c(3.1,4.1,0.6,0.6))
      plot(x=c(0,1),y=c(0,1.05*max(h$counts)),type="n",
           xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
      rect(xleft=par("usr")[1],xright=base$estimate,
           ybottom=par("usr")[3],ytop=par("usr")[4],
           border=NA,col=adjustcolor("dodgerblue3",alpha=0.3))
      rect(xleft=base$estimate,xright=par("usr")[2],
           ybottom=par("usr")[3],ytop=par("usr")[4],
           border=NA,col=adjustcolor("darkorange1",alpha=0.3))
      abline(v=seq(0,1,by=0.1),lty=3)
      #Plot colored histogram
      midpoints <- which(h$breaks==floor(100*base$estimate)/100 |
                          h$breaks==ceiling(100*base$estimate)/100)
      points.below <- which(h$breaks<=floor(100*base$estimate)/100)[-midpoints]
      rect(xleft=h$breaks[points.below],xright=h$breaks[points.below+1],
           ybottom=rep(0,times=length(points.below)),
           ytop=h$counts[points.below],
           col="dodgerblue3",border="dodgerblue4")
      points.above <- which(h$breaks>=ceiling(100*base$estimate)/100)[-midpoints]
      rect(xleft=h$breaks[points.above],xright=h$breaks[points.above+1],
           ybottom=rep(0,times=length(points.above)),
           ytop=h$counts[points.above],
           col="darkorange1",border="darkorange3")
      rect(xleft=h$breaks[min(midpoints)],xright=base$estimate,
           ybottom=0,ytop=h$counts[min(midpoints)],
           col="dodgerblue3",border="dodgerblue4")
      rect(xleft=base$estimate,xright=h$breaks[max(midpoints)],
           ybottom=0,ytop=h$counts[max(midpoints)],
           col="darkorange1",border="darkorange3")
      #Add marker & density curve
      dens <- density(res[,2])
      scaler <- max(h$counts)/max(dens$y)
      points(x=dens$x,type="l",y=scaler*dens$y)
      abline(v=base$estimate,lwd=2)
      #Add text
      text(x=c(par("usr")[1],par("usr")[2]),
           y=0.9*par("usr")[4],
           labels=c("Permutation Less Significant\nThan Observed Value",
                    "Permutation More Significant\nThan Observed Value"),
           pos=c(4,2),font=4,
           col=c("dodgerblue4","darkorange3"))
      if(base$estimate>=0.5){
        pos=2
      }else{
        pos=4
      }
      text(x=base$estimate,y=0.9*par("usr")[4],
           labels="Observed\nValue",pos=pos,font=4)
      #Add axes & labels
      axis(1,at=seq(0,1,0.2),labels=NA)
      axis(1,at=seq(0,1,0.2),labels=paste(seq(0,100,20),"%",sep=""),tick=F,line=-0.5)
      axis(2,las=2,cex.axis=0.75)
      mtext(1,text="Pct. of Bins with Element",line=1.5)
      mtext(2,text="Permutations (Count)",line=2.5)
      dev.off()
    }
  }

  #Print results to file
  if(measure=="mean"){
    results <- data.frame("mean.obs"=base$estimate[1],
                          "mean.perm"=est,
                          "SEM.perm"=std.err,
                          "p.perm"=p.emp)
  }else{
    results <- data.frame("pct.obs"=base$estimate[1],
                          "pct.perm"=est,
                          "SEM.perm"=std.err,
                          "p.perm"=p.emp)
  }
  write.table(results,paste(OUTDIR,"/",prefix,".annoBurden_results.txt",sep=""),
              col.names=T,row.names=F,quote=F,sep="\t")
}
