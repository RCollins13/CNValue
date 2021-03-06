#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate swarm plot of genic burden ORs by phenotype for Fig3a

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")
cols.ordered <- c(rep(cols.GERM[1],2),
                  rep(cols.NEURO[1],10),
                  rep(cols.SOMA[1],10),
                  rep(cols.CNCR[1],13))

#####Load data
pvals <- lapply(c("CNV","DEL","DUP"),function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",CNV,
                          "_E2_all_exonic.pvals.txt",sep=""),header=F)
  dat[,-1] <- apply(dat[,-1],2,as.numeric)
  colnames(dat) <- c("geneSet",phenos)
  return(dat)
})

#####Function to generate p-value scatter
plotPvals <- function(frame,row,CNV){
  #Get p-values
  vals <- unlist(-log10(pvals[[frame]][row,-1]))
  vals[which(is.infinite(vals))] <- 30

  #Prepare plotting area
  par(mar=c(1,4,1,0),bty="n")
  plot(x=1:35,y=vals,xlim=c(1,35),ylim=c(0,max(vals,na.rm=T)),
       pch=21,bg=cols.ordered,lwd=0.5,cex=1.35,
       xaxt="n",yaxt="n",xlab="",ylab="",
       panel.first=c(rect(xleft=c(par("usr")[1],2.5,12.5,22.5),
                          xright=c(2.5,12.5,22.5,par("usr")[2]),
                          ybottom=par("usr")[3],ytop=par("usr")[4],
                          border=NA,col=adjustcolor(c(cols.GERM[4],cols.NEURO[4],
                                                      cols.SOMA[4],cols.CNCR[4]),
                                                    alpha=0.4)),
                     abline(v=c(2.5,12.5,22.5),col=cols.CTRL[4]),
                     abline(h=unique(c(0,round(axTicks(2),0))),col=cols.CTRL[3]),
                     abline(h=-log10(0.05/nrow(pvals[[1]])),lty=2)))

  #Y-axis
  if(max(axTicks(2)>=4)){
    axis(2,at=unique(c(0,round(axTicks(2),0))),
         labels=NA)
    axis(2,at=unique(c(0,round(axTicks(2),0))),tick=F,
         labels=unique(c(0,round(axTicks(2),0))),
         line=-0.4,las=2)
  }else{
    axis(2,at=unique(c(0,round(axTicks(2),0),4)),
         labels=NA)
    axis(2,at=unique(c(0,round(axTicks(2),0),4)),tick=F,
         labels=unique(c(0,round(axTicks(2),0),4)),
         line=-0.4,las=2)
  }

  mtext(2,line=1.5,cex=0.7,
        text=bquote(paste(.(CNV),~~-log[10](italic(p)))))


  #Closing line
  axis(1,at=1:35,labels=NA,tck=0)

}



#Head/Neck in ID - CNV
plotPvals(1,108)

#3x2 plot - first column = specific, second column = shared
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/pheno_specific_geneSets.dotplots.pdf",sep=""),
    width=4,height=3.8)
par(mfrow=c(3,2))
#Axon guidance genes in ASD/Seizures - DEL
plotPvals(2,212,"DEL")
#Histone binding genes in ID, CRNL, CBST, and CLIV - CNV
plotPvals(1,178,"CNV")
#cartilage development in integument defects - CNV
plotPvals(1,240,"CNV")
#Transcriptional regulators in PSYCH/SCZ/BEHAV, CGEN, and CEND
plotPvals(3,172,"DUP")
#Highly expressed Endocrine genes in EMI, CNCR, CGEN, CEND, CLIV - DUP
plotPvals(3,93,"DUP")
#snoRNAs in NEURO, NDD, DD, SEIZ, HYPO, MSC, CBLD - CNV
plotPvals(1,161,"CNV")
dev.off()


#3x1 plot
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure3/pheno_specific_geneSets.dotplots.3_by_1.pdf",sep=""),
    width=3,height=3.8)
par(mfrow=c(3,1))
#cartilage development in integument defects - CNV
plotPvals(1,240,"CNV")
#ASD genes in ASD - DEL
plotPvals(2,21,"DEL")
#LCL genes in blood cancer - CNV
plotPvals(1,78,"CNV")
dev.off()


#Hematopoietic cell diff genes in liver cancer - DEL
plotPvals(2,239,"DEL")
#ASD genes in ASD - DEL
plotPvals(2,21,"DEL")
#Brain cancer associated genes - DUP
plotPvals(3,51,"DUP")







#Stem cell maintenance in 7 cancers
plotPvals(1,225,"CNV")

#Test
plotPvals(2,205,"DEL")

#Protein Phosphorylation in SEIZ HEAD GRO SKEL - DEL
plotPvals(2,205,"DEL")

#cartilage development in integument defects - CNV
plotPvals(1,240,"CNV")
#LCL genes in blood cancer - CNV
plotPvals(1,78,"CNV")
#Axon guidance genes in ASD/Seizures - DEL
plotPvals(2,212,"DEL")
#Ras Signalling in CSKN, CRNL, CBRN, CLNG, CHNK - DEL
plotPvals(2,257,"DEL")
#Nervous system dev genes in seizures - DUP
plotPvals(3,228,"DUP")
#Highly expressed Endocrine genes in EMI, CNCR, CGEN, CEND, CLIV - DUP
plotPvals(3,93,"DUP")

