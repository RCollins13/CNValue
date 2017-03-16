#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#fancyQQ: a dressed-up adaptation of Stephen Turner's qq() function from the qqman package
#Credit: S.D. Turner, qqman: an R package for visualizing GWAS results using Q-Q and Manhattan plots
#http://dx.doi.org/10.1101/005165

fancyQQ <- function (pvector,        #vector of p-values
                     nominal=0.05,   #threshold for nominal significance
                     adjusted=1E-8   #threshold for adjusted significance (e.g. genome-wide)
){
  if (!is.numeric(pvector))
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) &
                       is.finite(pvector) & pvector < 1 & pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  signif.cols <- rev(heat.colors(200*(-log10(adjusted)+1)))
  o.log10 <- ceiling(200*o)
  o.log10[which(o.log10>(200*(-log10(adjusted)+1)))] <- 200*(-log10(adjusted)+1)
  o.cols <- signif.cols[o.log10]
  o.cex.scale <- seq(0.75,2.25,by=0.05)
  o.cex.steps <- seq(0,max(o),by=max(o)/length(o.cex.scale))
  o.cex <- sapply(o,function(val){
    return(o.cex.scale[max(which(val>o.cex.steps))])
  })
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(x=e,y=o,pch=21,bg=o.cols,lwd=0.4,cex=o.cex,xaxt="n",yaxt="n",xlab="",ylab="",
       xaxs="i",yaxs="i",xlim=c(0,1.1*max(e)),ylim=c(0,1.1*max(o)),
       panel.first=c(rect(xleft=par("usr")[1],xright=par("usr")[2],
                          ybottom=par("usr")[3],ytop=-log10(nominal),
                          col=adjustcolor("cyan4",alpha=0.4),border=NA),
                     rect(xleft=par("usr")[1],xright=par("usr")[2],
                          ybottom=-log10(nominal),ytop=-log10(adjusted),
                          col=adjustcolor("cyan3",alpha=0.3),border=NA),
                     rect(xleft=par("usr")[1],xright=par("usr")[2],
                          ybottom=-log10(adjusted),ytop=par("usr")[4],
                          col=adjustcolor("cyan",alpha=0.2),border=NA),
                     abline(h=-log10(c(nominal,adjusted)),lty=2,col=c("cyan4","cyan3")),
                     abline(0,1,lwd=2),
                     text(x=c(0,0),y=-log10(c(nominal,adjusted)),cex=0.8,
                          labels=c("Nominal","Adjusted"),font=3,pos=4),
                     text(x=par("usr")[2],y=1.05*par("usr")[2],pos=2,labels="Obs ~ Exp",font=3,cex=0.8)))
  axis(1,at=seq(0,ceiling(max(e)),by=max(1,round(ceiling(max(e))/6,0))))
  axis(2,at=seq(0,ceiling(max(o)),by=max(1,round(ceiling(max(o))/6,0))),las=2)
  mtext(1,text=expression(Expected ~ ~-log[10](italic(p))),line=2.1)
  mtext(2,text=expression(Observed ~ ~-log[10](italic(p))),line=2.1)
}
