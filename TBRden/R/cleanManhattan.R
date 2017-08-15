#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#cleanManhattan: a dressed-up manhattan plot

cleanManhattan <- function (df,              #three-column data frame of chrom, pos, p-value
                            theme="green",   #select color theme (options: blue, green, red, orange, purple, brown)
                            nominal=0.05,    #threshold for nominal significance
                            adjusted=1E-8    #threshold for adjusted significance (e.g. genome-wide)
){
  #Set color theme vectors
  if(theme=="gray"){
    colors <- c("#6E6F70","#DCDDDF",NA,"black")
  }else if(theme=="purple"){
    colors <- c("#7B2AB3","#B07FD1",NA,"black")
  }else if(theme=="blue"){
    colors <- c("#00BFF4","#66D9F8",NA,"black")
  }else if(theme=="pink"){
    colors <- c("#EC008D","#F466BB",NA,"black")
  }else if(theme=="yellow"){
    colors <- c("#FFCB00","#FFE066",NA,"black")
  }
  #Gather list of unique contigs
  contigs <- unique(df[,1])
  contigs <- contigs[which(!(is.na(contigs)))]
  #Create index data frame
  indexes <- as.data.frame(t(sapply(contigs,function(chr){
    return(c(chr,0,max(df[which(df[,1]==chr),2])))
  })))
  indexes$sum <- cumsum(indexes[,3])
  indexes$bg <- rep(colors[1:2],ceiling(nrow(indexes)/2))[1:nrow(indexes)]
  #Create new plotting df with modified coordinates & color assignments
  df.plot <- as.data.frame(t(apply(df,1,function(row){
    return(c(row[1],
             as.numeric(row[2]+indexes[as.numeric(row[1]),4]-indexes[as.numeric(row[1]),3]),
             as.numeric(row[3]),
             indexes[as.numeric(row[1]),5]))
  })))
  df.plot[,2] <- as.numeric(as.character(df.plot[,2]))
  df.plot[,3] <- as.numeric(as.character(df.plot[,3]))
  #Prepare plotting window
  ymax <- -log10(df[which(!(is.na(df[,3]))),3])
  ymax <- max(max(ymax[which(!(is.infinite(ymax)))]),-log10(adjusted)+2)
  par(mar=c(2.1,3.1,0.6,0.6))
  plot(x=c(0,max(indexes[,4])),y=c(0,1.1*ymax),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
  rect(xleft=rep(par("usr")[1],3),xright=rep(par("usr")[2],3),
       ybottom=c(par("usr")[3],-log10(nominal),-log10(adjusted)),
       ytop=c(-log10(nominal),-log10(adjusted),par("usr")[4]),
       border=NA,col=c("gray60","gray80","gray95"))
  abline(v=indexes[,4],col="gray90")
  #Plots points & thresholds
  points(df.plot[,2],-log10(df.plot[,3]),
         cex=0.6,pch=21,lwd=0.01,bg=as.character(df.plot[,4]))
  abline(h=-log10(c(nominal,adjusted)),col=colors[3:4])
  #Adds axes & titles
  axis(1,at=indexes[,4],labels=NA)
  midpoints <- sapply(1:length(indexes[,2]),function(i){
    return(mean(c(indexes[i,4],indexes[i,4]-indexes[i,3])))
  })
  axis(1,at=midpoints,labels=indexes[,1],tick=F,line=-1.5,cex.axis=0.35)
  mtext(1,text="Chromosome",line=0.5)
  y.at <- seq(0,ceiling(par("usr")[4]),by=ceiling(par("usr")[4]/6))
  axis(2,at=y.at,labels=NA)
  axis(2,at=y.at,tick=F,line=-0.3,labels=y.at,cex.axis=0.75,las=2)
  mtext(2,text=expression(-log[10](italic(p))),line=1.5)
}
