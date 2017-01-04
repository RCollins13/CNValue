#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#binomialHelper: helper script to compute statistics and render plot
#from binomial_annotation_test.sh

binomialHelper <- function(obs,            #Number of observed hits
                           n,              #Number of loci tested
                           perm,           #Single-column text file containing hits for all matched
                           alt="greater",  #Alternative hypothesis: test bins are (greater|less) than random. Two-tailed test: two.sided
                           plot=T,         #Option to generate barplot on output
                           OUTDIR=NULL,    #Output directory for writing results
                           prefix="TBRden" #Prefix to be appended to results and plot (if optioned)
){
  #Load libraries
  require(plotrix)

  #Load data
  perm <- as.vector(read.table(perm,header=F)[,1])

  #Compute p-value binomial test
  p <- binom.test(obs,n,p=mean(perm/n),alternative=alt)$p.value

  #Compute odds ratio
  OR <- (obs/(n-obs))/(mean(perm/(n-perm)))

  #Plot (if optioned)
  if(plot==T){
    #Prep graphics output
    pdf(paste(OUTDIR,"/",prefix,".binomial_permutation_results.pdf",sep=""),
        height=3,width=6)
    layout(matrix(c(1,2),byrow=T,ncol=2),widths=c(4,1))
    #Prep plotting area for histogram
    xrange <- range(c(obs/n,perm/n))
    xrange <- c(max(c(xrange[1]-(diff(xrange)/3),0)),
                min(c(xrange[2]+(diff(xrange)/3),1)))
    breaks <- seq(xrange[1],xrange[2],by=diff(xrange)/50)
    h <- hist(perm/n,breaks=breaks,plot=F)
    par(mar=c(3.1,3.1,0.6,0.6))
    plot(x=xrange,y=c(0,1.1*max(h$counts)),type="n",
         xlab="",xaxt="n",xaxs="i",ylab="",yaxt="n",yaxs="i")
    rect(xleft=par("usr")[1],xright=par("usr")[2],
         ybottom=par("usr")[3],ytop=par("usr")[4],
         col="gray98")
    abline(v=seq(0,par("usr")[2],by=round(par("usr")[2]/6,3)),col="gray80")
    axis(1,at=seq(0,par("usr")[2],by=round(par("usr")[2]/6,3)),labels=NA)
    axis(1,at=seq(0,par("usr")[2],by=round(par("usr")[2]/6,3)),tick=F,line=-0.4,
         labels=paste(round(100*seq(0,par("usr")[2],by=round(par("usr")[2]/6,3)),1),"%",sep=""))
    mtext(1,text="Annotated Loci",line=1.8)
    abline(h=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),col="gray80")
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),labels=NA)
    axis(2,at=seq(0,par("usr")[4],by=round(par("usr")[4]/6,0)),tick=F,las=2,line=-0.4)
    mtext(2,text="Permutations",line=2.1)
    #Observed
    if(alt=="greater"){
      rect(xleft=obs/n,xright=par("usr")[2],
           ybottom=par("usr")[3],ytop=par("usr")[4],
           col=adjustcolor("darkorange1",alpha=0.7),border=NA)
    }else{
      rect(xleft=par("usr")[1],xright=obs/n,
           ybottom=par("usr")[3],ytop=par("usr")[4],
           col=adjustcolor("darkorange1",alpha=0.7),border=NA)
    }
    #Histogram
    rect(xleft=h$breaks[1:(length(h$breaks)-1)],
         xright=h$breaks[2:length(h$breaks)],
         ybottom=0,ytop=h$counts,
         col="dodgerblue3",border="dodgerblue4")
    #Annotations
    abline(v=c(obs/n,mean(perm/n)),lwd=2)
    points(x=c(obs/n,mean(perm/n)),y=rep(0.925*par("usr")[4],2),
           pch=23,cex=3,col="black",lwd=2,bg="white")
    text(x=c(obs/n,mean(perm/n)),y=rep(0.925*par("usr")[4],2),
         labels=c("O","E"),font=2,col=c("darkorange1","dodgerblue3"))
    #Scale results to permuted mean=1
    perm.scaled <- perm/mean(perm)
    obs.scaled <- obs/mean(perm)
    #Calculate 95% CI for perm.scaled
    perm.stde <- 2*std.error(perm.scaled)
    perm.CI <- c(mean(perm.scaled)-perm.stde,
                 mean(perm.scaled)+perm.stde)
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
    segments(x0=c(0.5,0.5,0.6),x1=c(0.7,0.7,0.6),
             y0=c(perm.CI[c(1,2,1)]),y1=c(perm.CI[c(1,2,2)]))
    #Category labels
    axis(1,at=0.6,tick=F,line=-1,labels=c("E"),font.axis=2,col.axis="dodgerblue3")
    axis(1,at=1.4,tick=F,line=-1,labels=c("O"),font.axis=2,col.axis="darkorange1")
    #Close device
    dev.off()
  }

  #Return results if OUTDIR isn't NULL
  if(!(is.null(OUTDIR))){
    res <- data.frame("Loci_Tested"=n,"Observed_Hits"=obs,"Expected_Hits"=mean(perm),
                      "Odds_Ratio"=OR,"p_Value"=p)
    write.table(res,paste(OUTDIR,"/",prefix,".TBRden_binomial_annotation_test.results.txt",sep=""),
                col.names=T,row.names=F,quote=F,sep="\t")
  }
}
