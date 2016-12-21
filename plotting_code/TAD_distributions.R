#!/usr/bin/env Rscript

#HST508 Final Project: Collins & Conway
#TAD intolerance & rare CNV maps

#Plotting code: TAD count by tissue per contig & TAD size by tissue




####Set parameters####
#Note: 21-color palette for tissue types adapted from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tissue.col=c("#771155","#AA4488","#CC99BB","#114477","#4477AA","#77AADD","#117777",
             "#44AAAA","#77CCCC","#117744","#44AA77","#88CCAA","#777711","#AAAA44",
             "#DDDD77","#774411","#AA7744","#DDAA77","#771122","#AA4455","#DD7788")
#Note: this vector does not include individual tissue replicates
tissues <- c("AD","AO","BL","CO","GM12878","H1","HC",
             "IMR90","LG","LI","LV","MES","MSC","NPC",
             "OV","PA","PO","RV","SB","SX","TRO")
options(scipen=6)
WRKDIR <- "/Users/collins/Desktop/RCollins/HMS/Courses/HST508/HST508_FinalProject/"
PLOTDIR <- paste(WRKDIR,"plots/",sep="")
#Create PLOTDIR if it doesn't already exist
if(!dir.exists(PLOTDIR)){
  dir.create(PLOTDIR)
}

####Load libraries####
require(vioplot)

####Plot tissue legends####
#Vertical (7x3)
pdf(paste(PLOTDIR,"tissue_legend.vertical.pdf",sep=""),height=2.75,width=3)
par(mar=rep(0.1,4))
plot(x=c(0,2.5),y=c(1,7),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
points(x=rep(0:2,times=7),
       y=sapply(1:7,rep,times=3),
       pch=19,col=tissue.col,cex=2)
text(x=rep(0:2,times=7),
     y=sapply(1:7,rep,times=3),
     pos=4,labels=tissues)
dev.off()
#Horizontal (3x7)
pdf(paste(PLOTDIR,"tissue_legend.horizontal.pdf",sep=""),height=1,width=7)
par(mar=rep(0.1,4))
plot(x=c(1,7.5),y=c(0.5,3.5),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
points(x=sapply(1:7,rep,times=3),
       y=rep(1:3,times=7),
       pch=19,col=tissue.col,cex=2)
text(x=sapply(1:7,rep,times=3),
     y=rep(1:3,times=7),
     pos=4,labels=tissues)
dev.off()

####Plot TAD count per contig by tissue####
#Read data
TAD.counts <- read.table(paste(WRKDIR,"plot_data/TAD_count_by_contig.txt",sep=""),header=T)
#Exclude replicates
TAD.counts <- TAD.counts[which(TAD.counts$tissue %in% tissues),]
#Prepare plot
pdf(paste(PLOTDIR,"TAD_counts_per_contig.pdf",sep=""),height=3,width=7)
par(mar=c(1.1,3.1,0.6,0.6))
#Plot
plot(x=c(0,23),y=c(0,225),type="n",ylab="",yaxt="n",xlab="",xaxt="n",yaxs="i",xaxs="i")
rect(xleft=seq(0,22,2),xright=seq(1,23,2),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray85",border=NA)
rect(xleft=seq(1,23,2),xright=seq(2,24,2),
     ybottom=par("usr")[3],ytop=par("usr")[4],
     col="gray95",border=NA)
abline(h=seq(0,200,50),col="gray60",lwd=0.5)
abline(h=seq(25,225,50),col="gray70",lwd=0.5,lty=2)
sapply(2:ncol(TAD.counts),function(i){
  points(x=seq(i-1.6,i-1.4,by=0.2/length(tissues))[1:21],
         y=TAD.counts[order(TAD.counts[,i]),i],
         pch=21,bg=tissue.col[order(TAD.counts[,i])],lwd=0.5,cex=0.7)
})
axis(1,at=0.5:22.5,labels=c(1:22,"X"),line=-2,tick=F,cex.axis=0.7,font=2)
axis(2,at=seq(0,200,50),las=2,labels=NA)
axis(2,at=seq(0,200,50),las=2,tick=F,labels=seq(0,200,50),line=-0.3,cex.axis=0.75)
mtext(1,text="Chromosome (GRCh37)",line=0,cex=0.75,font=2)
mtext(2,text="TADs per Contig (Count)",line=2,cex=0.75,font=2)
dev.off()

####Plot TAD size distribution by tissue####
#Read data
pdf(paste(PLOTDIR,"TAD_size_distributions.pdf",sep=""),height=3,width=7)
TAD.sizes <- sapply(tissues,function(tissue){
  data <- as.vector(read.table(paste(WRKDIR,"plot_data/TAD_size/",tissue,".TAD_size.txt",sep=""),header=F)[,1])
})
#Set log scaling
logvect <- log10(as.vector((sapply(4:8,function(x){return(1:9*(10^x))}))))
#Prep plot
par(mar=c(0.6,4.1,0.6,0.6))
plot(x=c(0,21),y=c(5,7),type="n",ylab="",yaxt="n",xlab="",xaxt="n",yaxs="i",xaxs="i")
abline(h=logvect,col="gray90",lwd=0.5)
abline(h=logvect[seq(5,39,9)],col="gray80")
abline(h=logvect[seq(1,39,9)],col="gray70",lwd=1.5)
axis(2,at=log10(c(100000,500000,1000000,5000000,10000000)),
     labels=c("100kb","500kb","1Mb","5Mb","10Mb"),las=2)
#Plot violins
sapply(1:length(TAD.sizes),function(i){
  vioplot(log10(TAD.sizes[[i]]),h=1/75,names=NA,col=tissue.col[i],at=i-0.5,add=T,pchMed=18)
})
dev.off()


