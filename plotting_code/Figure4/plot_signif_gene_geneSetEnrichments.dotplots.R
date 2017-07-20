#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate dotplots of rCNV signif gene overlaps with gene sets

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
phenos.GERM <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
                 "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
                 "MSC","EE","INT","EMI")
phenos.CNCR <- c("CNCR","CGEN","CSKN","CGST","CRNL",
                 "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")

#####Load libraries
require(plotrix)

############################################
#####Helper function to compute Fisher stats
############################################
fisherStats <- function(df){
  #Run Fisher tests
  fisher.res <- lapply(1:nrow(df),function(r){
    #Subset row
    vals <- as.numeric(df[r,-1])

    #CTRL
    CTRL.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                   vals[3]-vals[4],vals[4])),nrow=2)
    CTRL <- fisher.test(x=CTRL.df)
    #GERM
    GERM.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                   vals[5]-vals[6],vals[6])),nrow=2)
    GERM <- fisher.test(x=GERM.df)
    #NEURO
    NEURO.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                    vals[7]-vals[8],vals[8])),nrow=2)
    NEURO <- fisher.test(x=NEURO.df)
    #NDD
    NDD.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                  vals[9]-vals[10],vals[10])),nrow=2)
    NDD <- fisher.test(x=NDD.df)
    #PSYCH
    PSYCH.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                    vals[11]-vals[12],vals[12])),nrow=2)
    PSYCH <- fisher.test(x=PSYCH.df)
    #SOMA
    SOMA.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                   vals[13]-vals[14],vals[14])),nrow=2)
    SOMA <- fisher.test(x=SOMA.df)
    #CNCR
    CNCR.df <- matrix(as.numeric(c(vals[1]-vals[2],vals[2],
                                   vals[15]-vals[16],vals[16])),nrow=2)
    CNCR <- fisher.test(x=CNCR.df)

    #Collect output df
    res <- as.data.frame(t(data.frame("CTRL"=c(CTRL$estimate,CTRL$conf.int),
                                      "GERM"=c(GERM$estimate,GERM$conf.int),
                                      "NEURO"=c(NEURO$estimate,NEURO$conf.int),
                                      "NDD"=c(NDD$estimate,NDD$conf.int),
                                      "PSYCH"=c(PSYCH$estimate,PSYCH$conf.int),
                                      "SOMA"=c(SOMA$estimate,SOMA$conf.int),
                                      "CNCR"=c(CNCR$estimate,CNCR$conf.int))))
    colnames(res) <- c("OR","CI.lower","CI.upper")

    #Return values
    return(res)
  })

  #Assign names to results
  names(fisher.res) <- df[,1]

  #Return results
  return(fisher.res)
}

##############################################
#####Helper function to compute binomial stats
##############################################
binomStats <- function(df){
  #Run binomial tests
  binom.res <- lapply(1:nrow(df),function(r){
    #Subset row
    vals <- as.numeric(df[r,-1])

    #CTRL
    CTRL <- binom.test(x=vals[4],n=vals[3],p=vals[2]/vals[1])
    CTRL <- as.numeric(c(CTRL$estimate/as.numeric(CTRL$null.value),
                         CTRL$conf.int[1]/as.numeric(CTRL$null.value),
                         CTRL$conf.int[2]/as.numeric(CTRL$null.value),
                         CTRL$p.value))
    #GERM
    GERM <- binom.test(x=vals[6],n=vals[5],p=vals[2]/vals[1])
    GERM <- as.numeric(c(GERM$estimate/as.numeric(GERM$null.value),
                         GERM$conf.int[1]/as.numeric(GERM$null.value),
                         GERM$conf.int[2]/as.numeric(GERM$null.value),
                         GERM$p.value))
    #NEURO
    NEURO <- binom.test(x=vals[8],n=vals[7],p=vals[2]/vals[1])
    NEURO <- as.numeric(c(NEURO$estimate/as.numeric(NEURO$null.value),
                         NEURO$conf.int[1]/as.numeric(NEURO$null.value),
                         NEURO$conf.int[2]/as.numeric(NEURO$null.value),
                         NEURO$p.value))
    #NDD
    NDD <- binom.test(x=vals[10],n=vals[9],p=vals[2]/vals[1])
    NDD <- as.numeric(c(NDD$estimate/as.numeric(NDD$null.value),
                          NDD$conf.int[1]/as.numeric(NDD$null.value),
                          NDD$conf.int[2]/as.numeric(NDD$null.value),
                          NDD$p.value))
    #PSYCH
    PSYCH <- binom.test(x=vals[12],n=vals[11],p=vals[2]/vals[1])
    PSYCH <- as.numeric(c(PSYCH$estimate/as.numeric(PSYCH$null.value),
                        PSYCH$conf.int[1]/as.numeric(PSYCH$null.value),
                        PSYCH$conf.int[2]/as.numeric(PSYCH$null.value),
                        PSYCH$p.value))
    #SOMA
    SOMA <- binom.test(x=vals[14],n=vals[13],p=vals[2]/vals[1])
    SOMA <- as.numeric(c(SOMA$estimate/as.numeric(SOMA$null.value),
                          SOMA$conf.int[1]/as.numeric(SOMA$null.value),
                          SOMA$conf.int[2]/as.numeric(SOMA$null.value),
                          SOMA$p.value))
    #CNCR
    CNCR <- binom.test(x=vals[16],n=vals[15],p=vals[2]/vals[1])
    CNCR <- as.numeric(c(CNCR$estimate/as.numeric(CNCR$null.value),
                         CNCR$conf.int[1]/as.numeric(CNCR$null.value),
                         CNCR$conf.int[2]/as.numeric(CNCR$null.value),
                         CNCR$p.value))

    #Collect output df
    res <- as.data.frame(t(data.frame("CTRL"=CTRL,
                                      "GERM"=GERM,
                                      "NEURO"=NEURO,
                                      "NDD"=NDD,
                                      "PSYCH"=PSYCH,
                                      "SOMA"=SOMA,
                                      "CNCR"=CNCR)))
    colnames(res) <- c("fold","CI.lower","CI.upper","p")

    #Return values
    return(res)
  })

  #Assign names to results
  names(binom.res) <- df[,1]

  #Return results
  return(binom.res)
}

##############################################################
#####Helper function to read data & compute all relevant stats
##############################################################
readData <- function(CNV,VF,context,sig){
  #Read file
  df <- read.table(paste(WRKDIR,"plot_data/signif_genes_geneset_comparisons/subset/",
                         CNV,"_",VF,"_",context,"_",sig,".comparisons.txt",sep=""),
                   header=T,comment.char="")
  colnames(df)[1] <- "geneset"

  #Compute Fisher stats
  fisher.res <- fisherStats(df)

  #Compute binomial p-values
  binom.res <- binomStats(df)

  #Prepare output list
  results <- list("counts"=df,
                  "fisher"=fisher.res,
                  "binom"=binom.res)
  return(results)
}

########################################
#####Helper function to generate dotplot
########################################
plotDots <- function(dat,yaxis=T,
                     rows=NULL){
  #Subset data to columns of interest
  if(!is.null(rows)){
    dat$counts <- dat$counts[rows,]
    dat$binom <- dat$binom[rows]
    dat$fisher <- dat$fisher[rows]
  }

  #Prepare plot area
  par(mar=c(0.5,1,2.5,0.75),bty="n")
  plot(x=c(-0.1,5.1),y=c(0,-(nrow(dat$counts)-1)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i")

  #Draw gridlines
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=c(0:-nrow(dat$counts))-0.25,
       ytop=c(0:-nrow(dat$counts))+0.25,
       border=NA,col=cols.CTRL[4])
  abline(v=seq(0,4.2,0.5),col=cols.CTRL[3],lty=2)
  abline(v=seq(0,4.2,1),col=cols.CTRL[2],lty=2)
  abline(v=c(0,1))

  #Plot 95% CIs
  sapply(1:nrow(dat$counts),function(i){
    vals <- dat$binom[[i]]
    vals[which(vals[,1]==0),2:3] <- c(0,0)
    vals[which(vals[,2]>4.1),2] <- 4.1
    vals[which(vals[,3]>4.1),3] <- 4.1
    segments(x0=vals[c(1:3,6:7),2],
             x1=vals[c(1:3,6:7),3],
             y0=-(i-1)+seq(0.2,-0.2,-0.1),
             y1=-(i-1)+seq(0.2,-0.2,-0.1),
             lwd=0.5)
  })

  #Add X-axis
  axis(3,at=seq(0,4,0.5),col=cols.CTRL[1],tck=-0.01,labels=NA)
  axis(3,at=seq(0,4,1),tck=-0.02,labels=NA)
  # axis(3,at=seq(0,4,1),tick=F,line=-0.7,
  #      labels=c(0:3,">4"),las=2)

  # #Add Y-Axis (if optioned)
  # if(yaxis==T){
  #   sapply(1:nrow(dat$counts),function(i){
  #     axis(2,at=-i,labels=NA,tck=-0.02)
  #   })
  # }

  #Drop axis break
  # axis.break(3,breakpos=4.05,style="gap")
  rect(xleft=4,xright=4.3,ybottom=par("usr")[3],ytop=par("usr")[4],
       col="white",border=NA)
  abline(v=c(4,4.3),lty=c(2,1))

  #Plot points
  sapply(1:nrow(dat$counts),function(i){
    #Round data >2.5 to 2.55
    vals <- dat$binom[[i]][c(1:3,6:7),1]
    vals[which(vals>4)] <- 4
    points(x=vals,y=-(i-1)+seq(0.2,-0.2,-0.1),
           bg=c(cols.CTRL[1],cols.GERM[1],
                cols.NEURO[1],cols.SOMA[1],
                cols.CNCR[1]),pch=21)
  })

  #Prepare for p-value bars
  abline(v=seq(4.5,5.1,0.2),col=cols.CTRL[3])
  axis(3,at=c(4.3,5.1),labels=NA,tck=-0.02)

  #Plot -log10(p) as bars
  sapply(1:nrow(dat$counts),function(i){
    #Adjust p-vals
    pvals <- -log10(as.data.frame(dat$binom[i])[c(1:3,6:7),4])
    pvals[which(pvals>4)] <- 4
    pvals <- pvals/5

    #Plot bars
    rect(xleft=4.3,xright=4.3+pvals,
         ybottom=-(i-1)+seq(0.25,-0.25,-0.1)[2:6],
         ytop=-(i-1)+seq(0.25,-0.25,-0.1)[1:5],
         col=c(cols.CTRL[1],cols.GERM[1],
               cols.NEURO[1],cols.SOMA[1],
               cols.CNCR[1]),
         lwd=0.3)
  })

  #Add significance bars
  abline(v=4.3+((-log10(0.05))/5),col="red")
  axis(3,at=4.3+((-log10(0.05))/5),col="red",labels=NA,tck=-0.02,lwd=2)
}

###################
#####Generate plots
###################
#Create output directory
PLOTDIR <- paste(WRKDIR,"rCNV_map_paper/Figures/Figure4/geneSetEnrichment_dotplots/",sep="")
if(!dir.exists(PLOTDIR)){
  dir.create(PLOTDIR)
}
#Generate plots
sapply(c("CNV","DEL","DUP"),function(CNV){
  sapply(c("E2","E3","E4","N1"),function(VF){
    sapply(c("exonic","wholegene"),function(context){
      sapply(c("nominally","FDR","Bonferroni"),function(sig){
        #Read data
        dat <- readData(CNV,VF,context,sig)

        #Open device
        pdf(paste(PLOTDIR,CNV,"_",VF,"_",context,"_",sig,
                  ".geneSet_enrichments.pdf",sep=""),
            width=2.5,height=4)

        #Generate dotplots

        plotDots(dat,rows=c(1,3,9,12,14,15,20,17,16,22))

        #Close device
        dev.off()
      })
    })
  })
})
