#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#directTest: statistical modeling for TBRden_direct_burden_test.sh

directTest <- function(observed,         #Observed dCNVs for all elements
                       permuted,         #Matrix of permuted dCNV values (rows=permutations, columns=elements)
                       tail="upper",     #Specify "upper" or "lower"
                       Bonferroni=T,     #Perform Bonferroni correction on permuted p-values
                       QQ=T,             #Generate QQ plot
                       manhattan=T,      #Generate manhattan plot
                       manMaxY=50,       #Max y-value for manhattan plot
                       manColor="green", #Color for manhattan plot
                       return=F,         #Return results as data frame
                       OUTDIR=NULL,      #Output directory for writing results
                       prefix="TBRden",  #Prefix to be appended to results and QQ plot
                       gzip=T            #Gzip output
){
  #Sanity check input files
  if(!(file.exists(observed))){
    stop(paste("Input file ",observed," does not exist.",sep=""))
  }
  if(!(file.exists(permuted))){
    stop(paste("Input file ",permuted," does not exist.",sep=""))
  }

  #Set test tail
  if(tail=="upper"){
    lower.tail <- F
  }else if(tail=="lower"){
    lower=tail <- T
  }else{
    stop("Parameter 'tail' must be either 'upper' or 'lower'.")
  }

  #Import observed data as data frame
  OBS <- read.table(observed,header=F,comment.char="",sep="\t")
  if(ncol(OBS) == 4){
    colnames(OBS) <- c("chr","start","end","dCNV")
  }else if(ncol(OBS) == 5){
    colnames(OBS) <- c("chr","start","end","element_ID","dCNV")
  }else{
    stop(paste("Input file ",observed," does not contain 4 or 5 columns.",sep=""))
  }

  #Set threshold for genome-wide significance
  if(Bonferroni==T){
    p.thresh <- 0.05/nrow(OBS)
  }else{
    p.thresh <- 0.05
  }

  #Import permuted dCNVs
  PERM <- read.table(permuted,header=F,comment.char="",sep="\t")

  #Fit distribution to each element & compute Z-score & p-value
  fit.results <- as.data.frame(t(sapply(1:ncol(PERM),function(i){
    perm.mean <- mean(PERM[,i])
    perm.sd <- sd(PERM[,i])
    Z <- (OBS[i,5]-perm.mean)/perm.sd
    p <- pnorm(Z,mean=0,sd=1,lower.tail=lower.tail)
    return(c(OBS[i,5],perm.mean,perm.sd,Z,p))
  })))
  names(fit.results) <- c("dCNV.obs","dCNV.mean","dCNV.sd","Z","p")
  if(Bonferroni==T){
    fit.results$p.adj <- p.adjust(fit.results$p,method="bonferroni")
  }else{
    fit.results$p.adj <- fit.results$p
  }

  #Merge results with original bins
  res <- cbind(OBS[,-ncol(OBS)],fit.results)

  #Optional: generate QQ plot
  if(QQ==T){
    pdf(paste(OUTDIR,"/",prefix,".QQ.pdf",sep=""),height=6,width=6)
    cleanQQ(fit.results[,5],adjusted=p.thresh)
    dev.off()
  }

  #Clean output data frame and write to file
  write.table(res,paste(OUTDIR,"/",prefix,".TBRden_direct_test_results.bed",sep=""),
              col.names=T,row.names=F,sep="\t",quote=F)
  if(gzip==T){
    system(paste("gzip -f ",OUTDIR,"/",prefix,".TBRden_direct_test_results.bed",sep=""))
  }

  #Optional: generate manhattan plot
  if(manhattan==T){
    df <- res[,c(1,2,which(colnames(res)=="p"))]
    names(df) <- c("CHR","BP","P")
    df$CHR <- as.numeric(as.character(df$CHR))
    df[which(df$P<=10^-manMaxY),3] <- 10^-manMaxY
    pdf(paste(OUTDIR,"/",prefix,".manhattan.pdf",sep=""),height=4,width=8)
    cleanManhattan(df,adjusted=p.thresh,theme=manColor)
    dev.off()
  }

  #Optional: return results
  if(return==T){
    return(res)
  }

}
