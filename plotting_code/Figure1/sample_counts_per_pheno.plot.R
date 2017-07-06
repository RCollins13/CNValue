#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate barplot of patient counts per phenotype

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
logvect.all <- c(0,log10(as.vector(sapply(1:4,function(x){return((1:9)*10^x)}))))
logvect.mids <- log10(as.vector(sapply(1:5,function(x){return(c(5,10)*10^x)})))
logvect.mains <- log10(as.vector(sapply(1:5,function(x){return(10^x)})))

#####Prep data
#Read file
df <- read.table(paste(WRKDIR,"plot_data/figure1/sample_counts_by_group.txt",sep=""),
                 header=F,comment.char="")
names(df) <- c("group","color","line","count","description")
#Add column for log10 counts
df$logCount <- log10(df$count)

#####Modify NDD and PSYCH to have new color values
df[which(df$group=="NDD" | df$group=="PSYCH"),2] <- "#66d9f8"

#####Plot
#Initialize device
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/sample_counts_per_pheno.barplot.pdf",sep=""),
    height=5,width=2.5)
par(mar=c(0.5,0.5,2.5,0),bty="n")

#Prepare plot area
plot(y=c(0,-nrow(df)+1),x=c(2,log10(250000)),
     type="n",xaxt="n",yaxt="n",xaxs="i",xlab="",ylab="")
abline(v=logvect.all,col="gray92")
abline(v=logvect.mids,col="gray85")
abline(v=logvect.mains,col="gray70")

#X-axis
axis(3,logvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
axis(3,logvect.mids[1:8],labels=NA,tck=-0.02,col="gray20",lwd=1)
axis(3,logvect.mains[1:8],labels=NA,tck=-0.025,
  col="black",lwd=1.5)
axis(3,at=logvect.mids[1:8],line=-0.5,tick=F,las=2,cex.axis=0.7,
     labels=c(50,100,500,"1k","5k","10k","50k","100k"))
mtext(3,line=1.5,text="Subjects Matching Phenotype",cex=0.7)

#Plot bars
rect(xleft=0,xright=df$logCount,
     ybottom=seq(-1,-nrow(df))+0.15,
     ytop=seq(0,-nrow(df)+1),-0.15,
     col=df$color,border=df$line)

#Add text labels
sapply(1:nrow(df),function(i){
  if(df[i,1]!="SKIP"){
    text(x=df$logCount[i]-0.1,y=-i+0.5,pos=4,cex=0.6,
         labels=prettyNum(df$count[i],big.mark=","))
  }
})

#Final lines
abline(v=2,lwd=2)

#Close device
dev.off()
