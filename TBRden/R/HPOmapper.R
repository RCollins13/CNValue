#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#HPOmapper: function to assign HPO codes based on phenotype strings

HPOmapper <- function(phenotypes,      #Path to list of phenotypes
                      HPOmappings,     #Path to HPO map file. Column 1: keyword; Column 2: HPO codes
                      empty="0000000", #Automatically generate QQ plot
                      OUTFILE=NULL,    #Output file
                      gzip=T           #Gzip output file
){
  #Sanity check input files
  if(!(file.exists(phenotypes))){
    stop(paste("Input file ",phenotypes," does not exist.",sep=""))
  }
  if(!(file.exists(HPOmappings))){
    stop(paste("Input file ",HPOmappings," does not exist.",sep=""))
  }

  #Strings not as factors
  options(stringsAsFactors=F)

  #Import phenotypes
  phenos <- read.table(phenotypes,header=F,sep="\t")

  #Import HPO map
  HPOmap <- read.table(HPOmappings,header=F,sep="\t",quote="")

  #Consolidate HPO map
  HPOs <- unique(unlist(strsplit(HPOmap[,2],split=",")))
  HPOreg <- sapply(HPOs,function(HPO){
    terms <- unique(HPOmap[grep(HPO,as.character(HPOmap[,2])),1])
    terms <- gsub(":","",terms,fixed=T)
    terms <- gsub("/","",terms,fixed=T)
    terms <- gsub("\\","",terms,fixed=T)
    terms <- gsub("(","",terms,fixed=T)
    terms <- gsub(")","",terms,fixed=T)
    terms <- gsub("'","",terms,fixed=T)
    terms <- gsub("<","",terms,fixed=T)
    terms <- gsub(">","",terms,fixed=T)
    terms <- gsub(".","",terms,fixed=T)
    terms <- gsub("~","",terms,fixed=T)
    terms <- gsub("!","",terms,fixed=T)
    terms <- gsub("@","",terms,fixed=T)
    terms <- gsub("#","",terms,fixed=T)
    terms <- gsub("$","",terms,fixed=T)
    terms <- gsub("%","",terms,fixed=T)
    terms <- gsub("^","",terms,fixed=T)
    terms <- gsub("&","",terms,fixed=T)
    terms <- gsub("*","",terms,fixed=T)
    terms <- gsub("{","",terms,fixed=T)
    terms <- gsub("}","",terms,fixed=T)
    terms <- gsub("+","",terms,fixed=T)
    terms <- gsub("?","",terms,fixed=T)
    terms <- gsub(",","",terms,fixed=T)
    terms <- gsub(";","",terms,fixed=T)
    terms <- gsub("[","",terms,fixed=T)
    terms <- gsub("]","",terms,fixed=T)
    terms <- gsub("-","",terms,fixed=T)
    terms <- gsub("=","",terms,fixed=T)
    if(length(terms)>500){
      splits <- split(terms,ceiling(seq_along(terms)/500))
      regex <- lapply(splits,paste,collapse="|")
    }else{
      regex <- paste(terms,collapse="|")
    }
    return(regex)
  })

  #Make results df
  res <- data.frame("HPOs"=NA,"Phenos"=phenos[,1])

  #Iterate over unique HPO assignments
  for(i in 1:length(HPOs)){
    matches <- unique(unlist(lapply(HPOreg[[i]],function(regex){
      matches <- grep(as.character(regex),as.character(res[,2]),fixed=F,perl=T)
      return(matches)
    })))
    if(length(matches)>0){
      newHPOs <- sapply(matches,function(k){
        if(is.na(res[k,1])){
          return(HPOs[i])
        }else{
          return(paste(unique(c(unlist(strsplit(res[k,1],split=",")),HPOs[i])),collapse=","))
        }
      })
      res[matches,1] <- newHPOs
    }
  }

  #Assign empty HPO value to NA's
  res[which(is.na(res[,1])),1] <- empty

  #Write results to file
  if(is.null(OUTFILE)){
    return(res)
  }else{
    colnames(res)[1] <- "#HPOs"
    write.table(res,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
  }
  if(gzip==T){
    system(paste("gzip -f ",OUTFILE,sep=""))
  }

}
