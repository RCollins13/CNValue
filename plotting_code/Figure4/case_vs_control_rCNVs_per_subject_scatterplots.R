#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate case-vs-control rCNVs per gene scatterplots

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

#####Helper function to read data
readData <- function(pheno,VF,context){
  lapply(list("CNV","DEL","DUP"),function(CNV){
    #Read table
    dat <- read.table(paste(WRKDIR,"plot_data/figure4/geneScore_results/",
                            pheno,"_",CNV,"_",VF,"_",context,
                            ".geneScore_stats.txt",sep=""),
                      comment.char="",header=T)
    #Replace first column header
    names(dat)[1] <- "gene"
    #Return data
    return(dat)
  })
}

#####Helper function to plot individual scatterplot
CNVscatter <- function(df,color="red",maxval=NULL,
                       xaxis=T,yaxis=T,mar=c(3,3,0.5,0.5)){
  #Get max value between cases & controls
  if(is.null(maxval)){
    maxval <- max(c(df$control_ratio,df$case_ratio))
  }

  #Set plot parameters
  par(mar=mar)

  #Prepare plot area
  plot(df$control_ratio,df$case_ratio,type="n",
       xlim=c(0-(0.02*maxval),1.02*maxval),ylim=c(0-(0.02*maxval),1.02*maxval),
       xaxt="n",yaxt="n",xlab="",ylab="",
       xaxs="i",yaxs="i")

  #Shade nonsignificant middle
  cutoff <- qnorm(1-c(0.05/nrow(df)))/(mean(df$theta_Zscore)/mean(df$theta))
  polygon(x=c(par("usr")[1],cutoff-0.02*maxval,
              par("usr")[2],par("usr")[2],
              par("usr")[2]-cutoff,par("usr")[1]),
          y=c(par("usr")[3],par("usr")[3],
              par("usr")[4]-cutoff,par("usr")[4],
              par("usr")[4],cutoff-0.02*maxval),
          col=cols.CTRL[3],border=NA)

  #Add gridlines
  abline(h=axTicks(1),v=axTicks(2),col=cols.CTRL[4])

  #Add y~x line
  abline(0,1)

  #Add all points - black with low alpha
  points(df$control_ratio,df$case_ratio,cex=0.5,col=cols.CTRL[1])

  #Add lines indicating Bonferroni significance & shade nonsig middle
  abline(a=cutoff,b=1,lty=2)
  abline(a=-cutoff,b=1,lty=2)

  # #Add FDR significant points - red border, no fill
  # points(df$control_ratio[which(df$case_gt_control_ratio_FDR_q<=0.05)],
  #        df$case_ratio[which(df$case_gt_control_ratio_FDR_q<=0.05)],
  #        col="red",cex=0.5)
  # points(df$control_ratio[which(df$control_gt_case_ratio_FDR_q<=0.05)],
  #        df$case_ratio[which(df$control_gt_case_ratio_FDR_q<=0.05)],
  #        col="blue",cex=0.5)

  #Add Bonferroni significant points - red border, red fill
  points(df$control_ratio[which(df$case_gt_control_ratio_Bonferroni_p<=0.05)],
         df$case_ratio[which(df$case_gt_control_ratio_Bonferroni_p<=0.05)],
         col=color,cex=0.5,pch=19)
  points(df$control_ratio[which(df$control_gt_case_ratio_Bonferroni_p<=0.05)],
         df$case_ratio[which(df$control_gt_case_ratio_Bonferroni_p<=0.05)],
         col=color,cex=0.5,pch=19)

  # #Add FDR significant points from Fisher Test - red fill
  # points(df$control_ratio[which(df$case_gt_control_Fisher_FDR_q<=0.05)],
  #        df$case_ratio[which(df$case_gt_control_Fisher_FDR_q<=0.05)],
  #        col="red",cex=0.5)
  # points(df$control_ratio[which(df$control_gt_case_Fisher_FDR_q<=0.05)],
  #        df$case_ratio[which(df$control_gt_case_Fisher_FDR_q<=0.05)],
  #        col="blue",cex=0.5)

  #Add x-axis
  if(xaxis==T){
    axis(1,at=(axTicks(1)+(0.5*(axTicks(1)[2]-axTicks(1)[1])))[1:(length(axTicks(1))-1)],
         tck=-0.011,col="#dcdddf",label=NA)
    axis(1,at=axTicks(1),labels=NA,tck=-0.02)
    axis(1,at=axTicks(1),tick=F,line=-0.5,las=2,cex.axis=0.85,
         labels=paste(round(100*axTicks(1),2),"%",sep=""))
    # mtext(1,text="Adjusted rCNVs per Control Genome",line=1.1)
  }

  #Add y-axis
  if(yaxis==T){
    axis(2,at=(axTicks(2)+(0.5*(axTicks(2)[2]-axTicks(2)[1])))[1:(length(axTicks(2))-1)],
         tck=-0.01,col="#dcdddf",label=NA)
    axis(2,at=axTicks(2),labels=NA,tck=-0.02)
    axis(2,at=axTicks(2),tick=F,line=-0.5,las=2,cex.axis=0.85,
         labels=paste(round(100*axTicks(2),2),"%",sep=""))
    # mtext(2,text="Adjusted rCNVs per Case Genome",line=1.1)
  }

  #Add finishing border
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA,border="black")
}

#####Read data
# GERM.exonic.dfs <- readData("GERM","E4","exonic")
# NEURO.exonic.dfs <- readData("NEURO","E4","exonic")
# NDD.exonic.dfs <- readData("NDD","E4","exonic")
# PSYCH.exonic.dfs <- readData("PSYCH","E4","exonic")
# SOMA.exonic.dfs <- readData("SOMA","E4","exonic")
# CNCR.exonic.dfs <- readData("CNCR","E4","exonic")

#####Plot scatters
#GERM E4 exonic CNV
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_CNV.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.exonic.dfs[[1]],color="#333333",maxval=0.0015)
# dev.off()
# #GERM E4 exonic DEL
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DEL.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.exonic.dfs[[2]],color="#ed2126",yaxis=F,maxval=0.0015)
# dev.off()
# #GERM E4 exonic DUP
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_exonic_DUP.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.exonic.dfs[[3]],color="#3b53a4",yaxis=F,maxval=0.0015)
# dev.off()
#
# #####Get number of case/control significant genes for GERM E4 exonic CNVs
# sapply(1:3,function(i){
#   case <- length(which(GERM.exonic.dfs[[i]]$case_gt_control_ratio_Bonferroni_p<=0.05))
#   ctrl <- length(which(GERM.exonic.dfs[[i]]$control_gt_case_ratio_Bonferroni_p<=0.05))
#   return(c(case,ctrl))
# })



# #GERM E4 wholegene CNV
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_wholegene_CNV.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.wholegene.dfs[[1]])
# dev.off()
# #GERM E4 wholegene DEL
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_wholegene_DEL.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.wholegene.dfs[[2]])
# dev.off()
# #GERM E4 wholegene DUP
# png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/GERM_E4_wholegene_DUP.",
#           "case_control_scatter.png",sep=""),
#     width=3,height=3,units="in",res=1000)
# CNVscatter(GERM.wholegene.dfs[[3]])
# dev.off()


#SOMA exonic E4 DEL (for geneScore schematic)
# SOMA.exonic.dfs <- readData("SOMA","E4","exonic")
png(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/SOMA_E4_exonic_CNV.",
          "case_control_scatter.png",sep=""),
    width=2.5,height=2.5,units="in",res=400)
CNVscatter(SOMA.exonic.dfs[[3]],color="#FF6A09",
           maxval=0.001,xaxis=F,yaxis=F,
           mar=c(0.3,0.3,0.2,0.2))
dev.off()


####################################################
#####Theta histogram from SOMA vs CTRL E4 exonic DUP
####################################################
# SOMA.exonic.dfs <- readData("SOMA","E4","exonic")
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/theta_Z_hist.pdf",sep=""),
    width=2,height=2)
par(mar=c(0.3,0.3,0.2,0.2))
hist(SOMA.exonic.dfs[[3]]$theta_Zscore,
     breaks=200,xlim=c(-2,3),col=cols.CTRL[1],
     yaxs="i",yaxt="n",ylab="",main="",xlab="")
axis(2,at=seq(0,4000,1000),labels=NA)
dev.off()



