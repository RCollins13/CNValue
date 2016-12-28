#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2016 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Add columns for cross-disease and cross-CNV comparisions from master merged burden testing

#Set params
options(scipen=1000)
diseases <- c("DD","SCZ","DD_SCZ","CNCR")
CNVs <- c("CNV","DEL","DUP")
filts <- c("all","coding","noncoding")

#Read arguments
#first: input
#second: output
args <- commandArgs(trailingOnly=T)

#Read table
x <- read.table(args[1],comment.char="",header=T)

#Gather summary metrics
merged_results <- as.data.frame(t(sapply(1:nrow(x),function(i){
  vals <- x[i,-c(1:20)]
  names(vals) <- colnames(x)[-c(1:20)]

  #Report if bin was significant by ANY of CNV/DEL/DUP per disease
  ANY.p <- as.vector(sapply(diseases,function(disease){
    terms <- as.vector(sapply(CNVs,function(CNV){
      unlist(sapply(filts,function(filt){
        return(c(paste(disease,CNV,filt,"obs_p",sep="."),
                 paste(disease,CNV,filt,"perm_p",sep=".")))
      }))
    }))
    ANYp <- unlist(lapply(list(c(1,7,13),c(2,8,14),c(3,9,15),
                               c(4,10,16),c(5,11,17),c(6,12,18),
                               seq(1,17,2),seq(2,18,2)),
                          function(l){
      d <- vals[which(names(vals) %in% terms[l])]
      if(!(all(is.na(d)))){
        return(min(d,na.rm=T))
      }else{
        return(NA)
      }
    }))
    return(ANYp)
  }))

  #Report if bin was significant by ANY of CNV/DEL/DUP for ANY disease
  ANY.ANY.p <- as.vector(sapply(CNVs,function(CNV){
      terms <- as.vector(unlist(sapply(filts,function(filt){
        return(c(paste(CNV,filt,"obs_p",sep="."),
                 paste(CNV,filt,"perm_p",sep=".")))
      })))
      terms.idx <- as.vector(sapply(terms,function(term){
        grep(term,names(vals))
      }))
      ANY.ANY.p <- unlist(lapply(list(1:4,5:8,9:12,13:16,17:20,21:24,
                                      c(1:4,9:12,17:20),
                                      c(5:8,13:16,21:24)),
                                 function(l){
                            d <- vals[terms.idx[l]]
                            if(!(all(is.na(d)))){
                              return(min(d,na.rm=T))
                            }else{
                              return(NA)
                            }
                          }))
      return(ANY.ANY.p)
    }))
  final.p <- as.vector(sapply(filts,function(filt){
    terms <- c(paste(filt,"obs_p",sep="."),
               paste(filt,"perm_p",sep="."))
    terms.idx <- as.vector(sapply(terms,function(term){
      grep(term,names(vals))
    }))
    final.p <- unlist(lapply(list(1:12,13:24),function(l){
      d <- vals[terms.idx[l]]
      if(!(all(is.na(d)))){
        return(min(d,na.rm=T))
      }else{
        return(NA)
      }
    }))
    return(final.p)
  }))
  global_min.obs_p <- min(final.p[c(1,3,5)])
  if(!(all(is.na(final.p[c(2,4,6)])))){
    global_min.perm_p <- min(final.p[c(2,4,6)],na.rm=T)
  }else{
    global_min.perm_p <- NA
  }
  final.p <- c(final.p,global_min.obs_p,global_min.perm_p)

  #Concatenate results and print
  return(c(ANY.p,ANY.ANY.p,final.p))
})))
newcolnames <- as.vector(sapply(diseases,function(disease){
  sapply(c(filts,"ANY_FILTER"),function(filt){
    sapply(c("obs_p","perm_p"),function(measure){
      return(paste(disease,"ANY_CNV",filt,measure,sep="."))
    })
  })
}))
newcolnames <- c(newcolnames,
                 as.vector(sapply(c(CNVs,"ANY_CNV"),function(CNV){
                   sapply(c(filts,"ANY_FILTER"),function(filt){
                     sapply(c("obs_p","perm_p"),function(measure){
                       return(paste("ANY_DISEASE",CNV,filt,measure,sep="."))
                     })
                   })
})))
names(merged_results) <- newcolnames
results_out <- cbind(x,merged_results)
names(results_out)[1] <- "#chr"
write.table(results_out,args[2],col.names=T,row.names=F,quote=F,sep="\t")


