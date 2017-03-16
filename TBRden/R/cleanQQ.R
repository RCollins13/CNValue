#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#cleanQQ: an adaptation of Stephen Turner's qq() function from the qqman package
#Credit: S.D. Turner, qqman: an R package for visualizing GWAS results using Q-Q and Manhattan plots
#http://dx.doi.org/10.1101/005165

#Helper function to draw confidence interval
#Credit: http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
qqconf <- function(p.expected){
  n <- length(p.expected)
  mpts<-matrix(nrow=n*2, ncol=2)
  for(i in seq(from=1, to=n)) {
    mpts[i,1]<- -log10((i-.5)/n)
    mpts[i,2]<- -log10(qbeta(0.975, i, n-i))
    mpts[n*2+1-i,1]<- -log10((i-.5)/n)
    mpts[n*2+1-i,2]<- -log10(qbeta(0.025, i, n-i))
  }
  mpts <- as.data.frame(mpts)
  return(mpts)
}

#Main clean QQ function
cleanQQ <-  function(pvector,          #vector of observed p-values
                     color="red",      #color to shade significant points
                     nominal=0.05,     #threshold for nominal significance
                     adjusted=NULL,    #threshold for adjusted significance (NULL=Bonferroni)
                     xlab.omit=F,      #omit x-axis label
                     ylab.omit=F,      #omit y-axis label
                     print.stats=T     #print lambda and K-S p-value
){
  #Enxure p-value vector is numeric
  if (!is.numeric(pvector)){
    stop("Input must be numeric.")
  }

  #Clean pvector
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) &
                       is.finite(pvector) & pvector < 1 & pvector > 0]

  #Sort
  pvector <- pvector[order(pvector)]

  #Instantiate expected p quantiles
  expected <- ppoints(length(pvector))

  #Instantiate qq confidence intervals
  conf.int <- qqconf(expected)

  #Compute lambda and ks p-values
  ks.p <- ks.test(pvector,"punif")$p.value
  lambda <- dchisq(median(pvector),df=1)/dchisq(median(expected),df=1)

  #Convert to -log10 scale
  pvector <- -log10(pvector)
  expected <- -log10(expected)

  #Set adjusted significance to Bonferroni correction
  if(is.null(adjusted)){
    adjusted <- nominal/length(pvector)
  }

  #Get which points to color
  colors <- rep("gray20",length(pvector))
  colors[which(pvector>=-log10(adjusted))] <- color
  bg <- rep(NA,length(pvector))
  bg[which(pvector>=-log10(adjusted))] <- color

  #Plot
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(x=expected,y=pvector,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xaxs="i",yaxs="i",xlim=c(0,1.1*max(expected)),ylim=c(0,1.1*max(pvector)))
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border=NA,col="gray98")
  polygon(x=conf.int[,1],y=conf.int[,2],col="gray70",border=NA)
  abline(v=axTicks(1),col="gray92")
  abline(h=axTicks(2),col="gray92")
  abline(h=-log10(adjusted),col=color)
  abline(h=-log10(0.05),col="gray60",lty=2)
  points(x=expected,y=pvector,pch=21,col=colors,lwd=1.5,bg=bg)
  abline(0,1)
  axis(1,at=axTicks(1),labels=NA)
  axis(1,at=axTicks(1),tick=F,line=-0.3)
  if(xlab.omit==F){
    mtext(1,text=expression(Expected ~ ~-log[10](italic(p))),line=2.1,cex=0.7)
  }
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.3,las=2)
  if(ylab.omit==F){
    mtext(2,text=expression(Observed ~ ~-log[10](italic(p))),line=2.1,cex=0.7)
  }
  if(print.stats==T){
    text(par("usr")[1],0.88*par("usr")[4],pos=4,
         labels=(paste("K-S p=",sprintf("%.3f",ks.p),"\n",sep="")))
    text(par("usr")[1],0.85*par("usr")[4],pos=4,font=2,
         labels=bquote(lambda ~ .(paste("=",sprintf("%.2f",lambda),sep=""))))
  }
}
