#Test code to play with noncoding element clustering based on jaccard indexes

#Read data
x <- read.table("~/scratch/cleaned_noncoding_loci.jaccard_matrix.txt",header=T)
rownames(x) <- x[,1]
x <- x[,-1]
x <- apply(x,2,as.numeric)
rownames(x) <- colnames(x)

#Plot heatmap
png("~/scratch/heatmap.test.png",height=2000,width=2000,res=300)
heatmap(x,Rowv=NA,Colv=NA,col=colorRampPalette(c("white","black"))(100))
dev.off()

#Histogram of non-zero jaccard indexes < 1
x.v <- as.vector(as.matrix(x))
hist(x.v[which(x.v>0 & x.v<1)],breaks=20,col="black")

#Function to return number of unique clusters from NN clustering based on min jac cutoff
clusterElements <- function(x,minJac=0.2){
  #Write vector of all element IDs
  elements <- rownames(x)
  #Iterate over all elements and return list of matching IDs
  matches <- lapply(elements,function(element){
    return(c(names(which(x[which(colnames(x)==as.character(element)),]>minJac))))
  })
  #Iterate over resulting list and find any matches based on element IDs
  #Find A in B, repeat 10 times
  for(i in 1:10){
    matches <- lapply(matches,function(setA){
      hits <- lapply(matches,function(setB){
        return(any(setB %in% setA))
      })
      return(which(unlist(hits)))
    })
    # print(length(unique(unlist(lapply(matches,paste,collapse="_")))))
  }
  return(matches)
}

#Cluster with minJac = 0.2
matches <- clusterElements(x,minJac=0.2)
matches <- lapply(matches,function(elements){
  hits <- sapply(elements,function(element){
    return(rownames(x)[as.numeric(element)])
  })
  hits <- paste(hits,collapse=";")
  return(hits)
})
matches <- unique(unlist(matches))
write.table(matches,"~/scratch/clustered_elements_regBlocks.list",sep="\t",quote=F,col.names=F,row.names=F)
