#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate heatmap of patient sharing between germline phenotype groups

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
require(shape)
heatmap.palette <- colorRampPalette(c("#08176b","#138934","#ffff59"))(1001)

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Set color palettes
# require(RColorBrewer)
# cols.neuro <- colorRampPalette(c("white","#00BFF4"))(101)
# cols.soma <- colorRampPalette(c("white","#EC008D"))(101)
# cols.mixed <- colorRampPalette(c("white","#7B2AB3"))(101)

#####Read data
df <- read.table(paste(WRKDIR,"plot_data/rCNVdbSummaryFigure/germline_case_overlap.matrix.txt",sep=""))
phenos <- df[,1]
totals <- sapply(1:nrow(df),function(i){return(df[i,i+1])})
df <- df[,-1]

#####Compute jaccard index & scale to 1-1001 corresponding to 0-20%
scaled.df <- sapply(1:ncol(df),function(c){
  sapply(1:nrow(df),function(r){
    jaccard <- df[r,c]/(df[r,r]+df[c,c]-df[r,c])
    jaccard.rounded <- round(5000*jaccard,0)+1
    if(jaccard.rounded>1001){
      jaccard.rounded <- 1001
    }
    return(jaccard.rounded)
  })
})
#
# #####Transform data to pct overlap
# scaled.df <- as.data.frame(t(sapply(1:nrow(df),function(row){
#   total <- totals[row]
#   scaled.vals <- sapply(1:nrow(df),function(col){
#     return(df[row,col]/total)
#   })
#   scaled.vals <- round(1000*scaled.vals,0)
#   scaled.vals[which(scaled.vals>1000)] <- 1000
#   scaled.vals <- scaled.vals+1
#   return(scaled.vals)
# })))

#####Set shaded colors for final analysis
borcols <- c(rep(cols.NEURO[1],2),
             cols.NEURO[3],cols.NEURO[1],
             rep(cols.NEURO[3],6),
             cols.SOMA[1],
             rep(cols.SOMA[3],9))

#####Prepare plot
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/rCNVdbSummaryFigure/Germline_subject_phenotype_overlap.heatmap.pdf",sep=""),
    width=4,height=4)
par(mar=rep(0.2,4),bty="n")

#####Set up plotting area
plot(x=-1:nrow(df),y=-nrow(df):1,
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
# abline(h=-19.5:-0.5,v=0.5:19.5)
# Arrows(x0=par("usr")[1],x1=-0.22,y0=-19.5:-0.5,y1=-19.5:-0.5,
#        arr.type="triangle",arr.length=0.1,arr.width=0.15)
rect(xleft=(1:nrow(df))-0.9,xright=(1:nrow(df))-0.1,
     ybottom=0.1,ytop=0.5,col=borcols)
rect(xleft=-0.5,xright=-0.1,ybottom=-(1:nrow(df))+0.1,ytop=-((1:nrow(df))-0.9),
     col=borcols)

#####Plot heatmap
sapply(1:nrow(scaled.df),function(row){
  sapply(1:ncol(scaled.df),function(col){
    if(row==col){
      rect(xleft=col-1,xright=col,
           ybottom=-(row-1),ytop=-row,
           border=NA,col="#DCDDDF")
      points(x=col-0.5,y=-row+0.5,pch=18,cex=0.75,col="#6E6F70")
    }else{
        rect(xleft=col-1,xright=col,
             ybottom=-(row-1),ytop=-row,
             border=NA,col=heatmap.palette[scaled.df[row,col]])
    }
  })
})

#####Draw clean-up grid
# segments(x0=rep(0,ncol(df)+1),x1=c(1:ncol(df),ncol(df)),
#          y0=0:-(nrow(df)+1),y1=0:-(nrow(df)+1),
#          col="#F1F1F2")
# segments(y0=rep(-(ncol(df)),ncol(df)+1),y1=-c(0,0:ncol(df)),
#          x0=0:(nrow(df)+1),x1=0:(nrow(df)+1),
#          col="#F1F1F2")
segments(x0=rep(0,ncol(df)+1),x1=rep(ncol(df),ncol(df)),
         y0=0:-nrow(df),y1=0:-nrow(df),
         col="#F1F1F2")
segments(y0=rep(0,ncol(df)+1),y1=-rep(ncol(df),ncol(df)),
         x0=0:(nrow(df)+1),x1=0:(nrow(df)+1),
         col="#F1F1F2")
segments(x0=c(0,10),x1=c(20,10),
         y0=c(-10,0),y1=c(-10,-20),
         lty=2)
# rect(xleft=0,xright=nrow(df),ybottom=-nrow(df),ytop=0)

#####Close device
dev.off()

#####Plot scales/key
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/rCNVdbSummaryFigure/Germline_subject_phenotype_overlap.heatmap.key.pdf",sep=""),
    width=0.7,height=5)
par(mar=c(0.2,0.2,0.2,1),bty="n")
plot(x=c(0,1),y=c(0,-1001),
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
rect(xleft=0,xright=1,ybottom=-1001:-1,ytop=-1000:0,
     border=NA,col=rev(heatmap.palette[seq(1,1001,1)]))
#Draw borders
sapply(seq(1,1001,200),function(row){
  rect(xleft=0:2,xright=1:3,
       ybottom=-row,ytop=-(row-1),
       col=NA,border="#F1F1F2")
})
axis(4,at=seq(-1001,-1,200),labels=NA)
rect(xleft=0,xright=3,ybottom=par("usr")[3],ytop=par("usr")[4],col=NA)
#Close device
dev.off()




