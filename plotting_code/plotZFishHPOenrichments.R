#!/usr/bin/env R

#rCNV Map Project
#Winter 2018
#Talkowski Lab & Collaborators

#Copyright (c) 2018 Ryan Collins
#Distributed under terms of the MIT License

#Code to run hypergeometric tests for zfish HPO term enrichments

#####Set parameters & load libraries
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)


#####Read data
DEL <- read.table(paste(WRKDIR,"plot_data/zfish_HPO_collection/DEL_HPO.txt",sep=""),header=T)
DUP <- read.table(paste(WRKDIR,"plot_data/zfish_HPO_collection/DUP_HPO.txt",sep=""),header=T)


#####Helper functions
#Perform a hypergeometric test for a single gene
hyperGene <- function(all,obs){
  #Iterate over all endophenotypes and perform hypergeometric tests
  sapply(2:length(obs),function(i){
    p <- phyper(q=as.numeric(obs[i]),
                m=as.numeric(all[i]),
                n=as.numeric(all[1]-all[i]),
                k=as.numeric(obs[1]),
                lower.tail=F)
    return(p)
  })
}
#Perform hypergeometric tests for all genes
hyperAll <- function(df){
  ps <- t(sapply(2:nrow(df),function(r){
    if(df[r,2]>0){
      return(hyperGene(all=df[1,-1],obs=df[r,-1]))
    }else{
      return(rep(1,times=ncol(df)-2))
    }
  }))
  colnames(ps) <- colnames(df)[-c(1:2)]
  rownames(ps) <- df[-1,1]
  return(ps)
}
#Plot results of hyper matrix
plotHyperMatrix <- function(mat,title=NULL){
  #Sort matrix
  mat <- mat[order(rownames(mat)),rev(order(colnames(mat)))]

  #Round matrix values
  mat <- apply(mat,2,function(vals){
    vals[which(vals>0.1)] <- 1
    vals[which(vals<=0.1 & vals>0.01)] <- 2
    vals[which(vals<=0.01 & vals>0.001)] <- 3
    vals[which(vals<=0.001 & vals>0.0001)] <- 4
    vals[which(vals<0.0001)] <- 5
    return(vals)
  })
  mat <- as.data.frame(mat)

  #Prepare plot area
  par(mar=c(1,13,5,1),bty="n")
  plot(x=c(0,nrow(mat)),y=c(0,ncol(mat)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i")

  #Add gridlines & labels
  abline(h=0:ncol(mat),v=0:nrow(mat),col="gray80")
  axis(2,at=c(1:ncol(mat))-0.5,tick=F,line=-0.9,labels=colnames(mat),las=2,cex.axis=0.8)
  axis(3,at=c(1:nrow(mat))-0.5,tick=F,line=-0.9,labels=rownames(mat),las=2,cex.axis=0.8)
  mtext(3,line=4,text=title,font=2)

  #Plot rects
  rectCols <- c("gray95","gray85",rev(heat.colors(5))[3:5])
  sapply(1:nrow(mat),function(row){
    sapply(1:ncol(mat),function(col){
      val <- mat[row,col]
      rect(xleft=row-0.5-(0.1*val),xright=row-0.5+(0.1*val),
           ybottom=col-0.5-(0.1*val),ytop=col-0.5+(0.1*val),
           border=NA,col=rectCols[val])
    })
  })
}


#####Run tests on data
DEL.mat <- hyperAll(DEL)
DUP.mat <- hyperAll(DUP)
pdf(paste(WRKDIR,"SSC_LCL_nanostring_and_Zebrafish_dupValidation/HPO_enrichments.pdf",sep=""),
    height=10,width=8)
par(mfrow=c(2,1))
plotHyperMatrix(DEL.mat,title="Deletion HPO Enrichments")
plotHyperMatrix(DUP.mat,title="Duplication HPO Enrichments")
dev.off()


#####Plot key
par(mar=c(0.1,0.1,0.1,6))
plot(x=c(0,1),y=c(0,5),type="n",
     xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
rectCols <- c("gray95","gray85",rev(heat.colors(5))[3:5])
rect(xleft=0.5-seq(0.1,0.5,0.1),
     xright=0.5+seq(0.1,0.5,0.1),
     ybottom=seq(0.5,4.5)-seq(0.1,0.5,0.1),
     ytop=seq(0.5,4.5)+seq(0.1,0.5,0.1),
     border="black",col=rectCols)
axis(4,at=0.5:4.5,tick=F,las=2,line=-0.8,
     labels=c(">0.1","0.1-0.01","0.01-0.001",
              "0.001-0.0001","<0.0001"))



