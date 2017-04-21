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

#####Set color palettes
require(RColorBrewer)
cols.neuro <- colorRampPalette(c("white","#00BFF4"))(101)
cols.soma <- colorRampPalette(c("white","#EC008D"))(101)
cols.mixed <- colorRampPalette(c("white","#7B2AB3"))(101)

#####Read data
df <- read.table(paste(WRKDIR,"plot_data/figure1/germline_case_overlap.matrix.txt",sep=""))
phenos <- df[,1]
totals <- sapply(1:nrow(df),function(i){return(df[i,i+1])})
df <- df[,-1]

#####Transform data to pct overlap
scaled.df <- as.data.frame(t(sapply(1:nrow(df),function(row){
  total <- totals[row]
  scaled.vals <- sapply(1:nrow(df),function(col){
    return(df[row,col]/total)
  })
  scaled.vals <- round(100*scaled.vals,0)
  scaled.vals[which(scaled.vals>100)] <- 100
  scaled.vals <- scaled.vals+1
  return(scaled.vals)
})))

#####Prepare plot
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/Germline_subject_phenotype_overlap.heatmap.pdf",sep=""),
    width=4,height=4)
par(mar=rep(0.2,4),bty="n")

#####Set up plotting area
plot(x=-1:nrow(df),y=-nrow(df):1,
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
# abline(h=-19.5:-0.5,v=0.5:19.5)
Arrows(x0=par("usr")[1],x1=-0.22,y0=-19.5:-0.5,y1=-19.5:-0.5,
       arr.type="triangle",arr.length=0.1,arr.width=0.15)
rect(xleft=(1:nrow(df))-0.9,xright=(1:nrow(df))-0.1,
     ybottom=0.1,ytop=0.5,col=c(rep("#00BFF4",10),rep("#EC008D",10)))
rect(xleft=-1,xright=-0.6,ybottom=-(1:nrow(df))+0.1,ytop=-((1:nrow(df))-0.9),
     col=c(rep("#00BFF4",10),rep("#EC008D",10)))

#####Plot heatmap
#NEURO
sapply(1:10,function(row){
  sapply(1:ncol(scaled.df),function(col){
    if(row==col){
      rect(xleft=col-1,xright=col,
           ybottom=-(row-1),ytop=-row,
           border=NA,col="#DCDDDF")
      points(x=col-0.5,y=-row+0.5,pch=18,cex=0.75,col="#6E6F70")
    }else{
      if(col>10){
        rect(xleft=col-1,xright=col,
             ybottom=-(row-1),ytop=-row,
             border=NA,col=cols.mixed[scaled.df[row,col]])
      }else{
        rect(xleft=col-1,xright=col,
             ybottom=-(row-1),ytop=-row,
             border=NA,col=cols.neuro[scaled.df[row,col]])
      }
    }
  })
})
#SOMA
sapply(11:nrow(scaled.df),function(row){
  sapply(1:ncol(scaled.df),function(col){
    if(row==col){
      rect(xleft=col-1,xright=col,
           ybottom=-(row-1),ytop=-row,
           border=NA,col="#DCDDDF")
      points(x=col-0.5,y=-row+0.5,pch=18,cex=0.75,col="#6E6F70")
    }else{
      if(col<11){
        rect(xleft=col-1,xright=col,
             ybottom=-(row-1),ytop=-row,
             border=NA,col=cols.mixed[scaled.df[row,col]])
      }else{
        rect(xleft=col-1,xright=col,
             ybottom=-(row-1),ytop=-row,
             border=NA,col=cols.soma[scaled.df[row,col]])
      }
    }
  })
})

#####Draw clean-up grid
sapply(1:nrow(df),function(i){
  rect(xleft=0,xright=nrow(df),
       ybottom=-i,ytop=-(i-1),
       col=NA,border="#F1F1F2")
  rect(xleft=i-1,xright=i,ybottom=0,ytop=-nrow(df),
       col=NA,border="#F1F1F2")
})
segments(x0=c(0,10),x1=c(20,10),
         y0=c(-10,0),y1=c(-10,-20))
# rect(xleft=0,xright=nrow(df),ybottom=-nrow(df),ytop=0)

#####Close device
dev.off()

#####Plot scales/key
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/Germline_subject_phenotype_overlap.heatmap.key.pdf",sep=""),
    width=0.7,height=5)
par(mar=c(0.2,0.2,0.2,1),bty="n")
plot(x=c(0,3),y=c(0,-5),
     type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
#NEURO
rect(xleft=0,xright=1,ybottom=-5:-1,ytop=-4:0,
     border=NA,col=cols.neuro[seq(1,101,20)])
#SOMA
rect(xleft=1,xright=2,ybottom=-5:-1,ytop=-4:0,
     border=NA,col=cols.soma[seq(1,101,20)])
#GERM
rect(xleft=2,xright=3,ybottom=-5:-1,ytop=-4:0,
     border=NA,col=cols.mixed[seq(1,101,20)])
#Draw borders
sapply(1:5,function(row){
  rect(xleft=0:2,xright=1:3,
       ybottom=-row,ytop=-(row-1),
       col=NA,border="#F1F1F2")
})
axis(4,at=-5:0,labels=NA)
rect(xleft=0,xright=3,ybottom=-5,ytop=0,col=NA)
#Close device
dev.off()




