#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#tHelper: helper script to compute t-test statistics and render plot
#from hotspot_annotation_test.sh

tHelper <- function(obs,                            #Single-column text file containing hits for all observed loci
                    perm,                           #Single-column text file containing mean number of hits per locus for all matched permutations
                    alt="greater",                  #Alternative hypothesis: test bins are (greater|less) than random. Two-tailed test: two.sided
                    plot=T,                         #Option to generate barplot on output
                    OUTDIR=NULL,                    #Output directory for writing results
                    prefix="TBRden_annotation_test" #Prefix to be appended to results and plot (if optioned)
){
  #Load libraries
  require(plotrix)

  #Load data
  obs <- as.vector(read.table(obs,header=F)[,1])
  perm <- as.vector(read.table(perm,header=F)[,1])

  #Compute p-value from t-test
  p <- t.test(obs,mu=mean(perm),alternative=alt)$p.value

  #Compute empirical permuted p-value
  if(alt=="greater"){
    p.perm <- (length(which(perm>=mean(obs)))+1)/(length(perm)+1)
  }else{
    p.perm <- (length(which(perm<=mean(obs)))+1)/(length(perm)+1)
  }

  #Compute fold-change
  fold <- mean(obs)/mean(perm)

  #Plot (if optioned)
  if(plot==T){
    #Prep graphics output
    pdf(paste(OUTDIR,"/",prefix,".ttest_permutation_results.pdf",sep=""),
        height=3,width=6)
    layout(matrix(c(1,2),byrow=T,ncol=2),widths=c(4,1))
    #Prep plotting area for histogram
    xrange <- range(c(mean(obs),perm))
    xrange <- c(max(c(xrange[1]-(diff(xrange)/3),0)),
                xrange[2]+(diff(xrange)/3))
    breaks <- seq(xrange[1],xrange[2],by=diff(xrange)/50)
    h <- hist(perm,breaks=breaks,plot=F)
    par(mar=c(3.1,3.1,0.6,0.6))
    plot(x=xrange,y=c(0,1.1*max(h$counts)),type="n",
         xlab="",xaxt="n",xaxs="i",ylab="",yaxt="n",yaxs="i")
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=par("usr")[3],ytop=par("usr")[4],
         col="gray98")
    abline(v=seq(0,par("usr")[2],by=round(par("usr")[2]/6,2)),col="gray80")
    axis(1,at=seq(0,par("usr")[2],by=round(par("usr")[2]/6,2)),labels=NA)
    axis(1,at=seq(0,par("usr")[2],by=round(par("usr")[2]/6,2)),tick=F,line=-0.4,
         labels=seq(0,par("usr")[2],by=round(par("usr")[2]/6,2)))
    mtext(1,text="Mean Elements per Locus",line=1.8)
    abline(h=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),col="gray80")
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),labels=NA)
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),tick=F,las=2,line=-0.4)
    mtext(2,text="Permutations",line=2.1)
    #Observed
    if(alt=="greater"){
      rect(xleft=mean(obs),xright=par("usr")[2],
           ybottom=par("usr")[3],ytop=par("usr")[4],
           col=adjustcolor("darkorange1",alpha=0.7),border=NA)
    }else{
      rect(xleft=par("usr")[1],xright=mean(obs),
           ybottom=par("usr")[3],ytop=par("usr")[4],
           col=adjustcolor("darkorange1",alpha=0.7),border=NA)
    }
    #Histogram
    rect(xleft=h$breaks[1:(length(h$breaks)-1)],
         xright=h$breaks[2:length(h$breaks)],
         ybottom=0,ytop=h$counts,
         col="dodgerblue3",border="dodgerblue4")
    #Annotations
    abline(v=c(mean(obs),mean(perm)),lwd=2)
    points(x=c(mean(obs),mean(perm)),y=rep(0.925*par("usr")[4],2),
           pch=23,cex=3,col="black",lwd=2,bg="white")
    text(x=c(mean(obs),mean(perm)),y=rep(0.925*par("usr")[4],2),
         labels=c("O","E"),font=2,col=c("darkorange1","dodgerblue3"))
    #Scale results to permuted mean=1
    perm.scaled <- perm/mean(perm)
    obs.scaled <- mean(obs)/mean(perm)
    #Calculate 95% CI for perm.scaled
    perm.sd <- 2*sd(perm.scaled)
    obs.CI <- c(obs.scaled-perm.sd,
                obs.scaled+perm.sd)
    #Prep plot area for barplot
    plot(x=c(0,2),y=c(0,1.1*max(c(perm.scaled,obs.scaled))),type="n",
         xlab="",xaxt="n",xaxs="i",ylab="",yaxt="n",yaxs="i")
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=par("usr")[3],ytop=par("usr")[4],
         col="gray98")
    abline(h=seq(0,par("usr")[4],by=round(par("usr")[4]/6,1)),col="gray80")
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,1)),labels=NA)
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,1)),tick=F,las=2,line=-0.4)
    abline(h=1,lwd=2)
    mtext(2,text="Fold-Change vs Expected",line=2.1)
    #Bars
    rect(xleft=c(0.2,1),xright=c(1,1.8),
         ybottom=c(0,0),ytop=c(1,obs.scaled),
         col=c("dodgerblue3","darkorange1"),
         border=c("dodgerblue4","darkorange2"))
    #95% CI
    segments(x0=c(1.3,1.3,1.4),x1=c(1.5,1.5,1.4),
             y0=c(obs.CI[c(1,2,1)]),y1=c(obs.CI[c(1,2,2)]))
    #Category labels
    axis(1,at=0.6,tick=F,line=-1,labels=c("E"),font.axis=2,col.axis="dodgerblue3")
    axis(1,at=1.4,tick=F,line=-1,labels=c("O"),font.axis=2,col.axis="darkorange1")
    #Close device
    dev.off()
  }

  #Return results if OUTDIR isn't NULL
  if(!(is.null(OUTDIR))){
    res <- data.frame("Loci_Tested"=n,"Observed_Mean"=mean(obs),"Expected_Mean"=mean(perm),
                      "Fold_Change"=fold,"Fold_Change_95pct_Lower"=obs.CI[1],"Fold_Change_95pct_Upper"=obs.CI[2],
                      "tTest_p"=p,"Perm_p"=p.perm)
    write.table(res,paste(OUTDIR,"/",prefix,".TBRden__annotation_test.results.txt",sep=""),
                col.names=T,row.names=F,quote=F,sep="\t")
  }
}
