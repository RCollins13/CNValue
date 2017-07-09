#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Exploratory code for per-gene CNV burden scores

#####Set parameters
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")
#Sample sizes (CTRL, GERM, CNCR)
nsamp <- c(38628,63629,10844)

#####Load required libraries
require(MASS)

#####Set color vectors
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","#B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")



#####Read data
#GERM
GERM.E2.exonic <- readData("GERM","E2","exonic")
GERM.E2.wholegene <- readData("GERM","E2","wholegene")
GERM.E3.exonic <- readData("GERM","E3","exonic")
GERM.E3.wholegene <- readData("GERM","E3","wholegene")
GERM.E4.exonic <- readData("GERM","E4","exonic")
GERM.E4.wholegene <- readData("GERM","E4","wholegene")
GERM.N1.exonic <- readData("GERM","N1","exonic")
GERM.N1.wholegene <- readData("GERM","N1","wholegene")
#NDD
NDD.E2.exonic <- readData("NDD","E2","exonic")
NDD.E2.wholegene <- readData("NDD","E2","wholegene")
NDD.E3.exonic <- readData("NDD","E3","exonic")
NDD.E3.wholegene <- readData("NDD","E3","wholegene")
NDD.E4.exonic <- readData("NDD","E4","exonic")
NDD.E4.wholegene <- readData("NDD","E4","wholegene")
NDD.N1.exonic <- readData("NDD","N1","exonic")
NDD.N1.wholegene <- readData("NDD","N1","wholegene")
#PSYCH
PSYCH.E2.exonic <- readData("PSYCH","E2","exonic")
PSYCH.E2.wholegene <- readData("PSYCH","E2","wholegene")
PSYCH.E3.exonic <- readData("PSYCH","E3","exonic")
PSYCH.E3.wholegene <- readData("PSYCH","E3","wholegene")
PSYCH.E4.exonic <- readData("PSYCH","E4","exonic")
PSYCH.E4.wholegene <- readData("PSYCH","E4","wholegene")
PSYCH.N1.exonic <- readData("PSYCH","N1","exonic")
PSYCH.N1.wholegene <- readData("PSYCH","N1","wholegene")
#ASD
ASD.E2.exonic <- readData("ASD","E2","exonic")
ASD.E2.wholegene <- readData("ASD","E2","wholegene")
ASD.E3.exonic <- readData("ASD","E3","exonic")
ASD.E3.wholegene <- readData("ASD","E3","wholegene")
ASD.E4.exonic <- readData("ASD","E4","exonic")
ASD.E4.wholegene <- readData("ASD","E4","wholegene")
ASD.N1.exonic <- readData("ASD","N1","exonic")
ASD.N1.wholegene <- readData("ASD","N1","wholegene")
#SEIZ
SEIZ.E2.exonic <- readData("SEIZ","E2","exonic")
SEIZ.E2.wholegene <- readData("SEIZ","E2","wholegene")
SEIZ.E3.exonic <- readData("SEIZ","E3","exonic")
SEIZ.E3.wholegene <- readData("SEIZ","E3","wholegene")
SEIZ.E4.exonic <- readData("SEIZ","E4","exonic")
SEIZ.E4.wholegene <- readData("SEIZ","E4","wholegene")
SEIZ.N1.exonic <- readData("SEIZ","N1","exonic")
SEIZ.N1.wholegene <- readData("SEIZ","N1","wholegene")

#####Code to plot CNV count correlations
#exonic DEL vs gene length
plot(GERM.N1.exonic[[2]]$gene_length,GERM.N1.exonic[[2]]$all_CNV)
plot(GERM.N1.exonic[[2]]$gene_length,GERM.N1.exonic[[2]]$case_CNV)
plot(GERM.N1.exonic[[2]]$gene_length,GERM.N1.exonic[[2]]$control_CNV)
#exonic DUP vs gene length
plot(GERM.N1.exonic[[3]]$gene_length,GERM.N1.exonic[[3]]$all_CNV)
plot(GERM.N1.exonic[[3]]$gene_length,GERM.N1.exonic[[3]]$case_CNV)
plot(GERM.N1.exonic[[3]]$gene_length,GERM.N1.exonic[[3]]$control_CNV)
#whole-gene DEL vs gene length
plot(GERM.N1.wholegene[[2]]$gene_length,GERM.N1.wholegene[[2]]$all_CNV)
plot(GERM.N1.wholegene[[2]]$gene_length,GERM.N1.wholegene[[2]]$case_CNV)
plot(GERM.N1.wholegene[[2]]$gene_length,GERM.N1.wholegene[[2]]$control_CNV)
#whole-gene DUP vs gene length
plot(GERM.N1.wholegene[[3]]$gene_length,GERM.N1.wholegene[[3]]$all_CNV)
plot(GERM.N1.wholegene[[3]]$gene_length,GERM.N1.wholegene[[3]]$case_CNV)
plot(GERM.N1.wholegene[[3]]$gene_length,GERM.N1.wholegene[[3]]$control_CNV)
#exonic DEL vs exonic bases
plot(GERM.N1.exonic[[2]]$exonic_bases,GERM.N1.exonic[[2]]$all_CNV,xlim=c(0,40000))
plot(GERM.N1.exonic[[2]]$exonic_bases,GERM.N1.exonic[[2]]$case_CNV,xlim=c(0,40000))
plot(GERM.N1.exonic[[2]]$exonic_bases,GERM.N1.exonic[[2]]$control_CNV,xlim=c(0,40000))
#exonic DUP vs exonic bases
plot(GERM.N1.exonic[[3]]$exonic_bases,GERM.N1.exonic[[3]]$all_CNV,xlim=c(0,40000))
plot(GERM.N1.exonic[[3]]$exonic_bases,GERM.N1.exonic[[3]]$case_CNV,xlim=c(0,40000))
plot(GERM.N1.exonic[[3]]$exonic_bases,GERM.N1.exonic[[3]]$control_CNV,xlim=c(0,40000))
#whole-gene DEL vs exonic bases
plot(GERM.N1.wholegene[[2]]$exonic_bases,GERM.N1.wholegene[[2]]$all_CNV,xlim=c(0,40000))
plot(GERM.N1.wholegene[[2]]$exonic_bases,GERM.N1.wholegene[[2]]$case_CNV,xlim=c(0,40000))
plot(GERM.N1.wholegene[[2]]$exonic_bases,GERM.N1.wholegene[[2]]$control_CNV,xlim=c(0,40000))
#whole-gene DUP vs exonic bases
plot(GERM.N1.wholegene[[3]]$exonic_bases,GERM.N1.wholegene[[3]]$all_CNV,xlim=c(0,40000))
plot(GERM.N1.wholegene[[3]]$exonic_bases,GERM.N1.wholegene[[3]]$case_CNV,xlim=c(0,40000))
plot(GERM.N1.wholegene[[3]]$exonic_bases,GERM.N1.wholegene[[3]]$control_CNV,xlim=c(0,40000))
#exonic DEL vs GC
plot(GERM.N1.exonic[[2]]$GC,GERM.N1.exonic[[2]]$all_CNV)
plot(GERM.N1.exonic[[2]]$GC,GERM.N1.exonic[[2]]$case_CNV)
plot(GERM.N1.exonic[[2]]$GC,GERM.N1.exonic[[2]]$control_CNV)
#exonic DUP vs GC
plot(GERM.N1.exonic[[3]]$GC,GERM.N1.exonic[[3]]$all_CNV)
plot(GERM.N1.exonic[[3]]$GC,GERM.N1.exonic[[3]]$case_CNV)
plot(GERM.N1.exonic[[3]]$GC,GERM.N1.exonic[[3]]$control_CNV)
#whole-gene DEL vs GC
plot(GERM.N1.wholegene[[2]]$GC,GERM.N1.wholegene[[2]]$all_CNV)
plot(GERM.N1.wholegene[[2]]$GC,GERM.N1.wholegene[[2]]$case_CNV)
plot(GERM.N1.wholegene[[2]]$GC,GERM.N1.wholegene[[2]]$control_CNV)
#whole-gene DUP vs GC
plot(GERM.N1.wholegene[[3]]$GC,GERM.N1.wholegene[[3]]$all_CNV)
plot(GERM.N1.wholegene[[3]]$GC,GERM.N1.wholegene[[3]]$case_CNV)
plot(GERM.N1.wholegene[[3]]$GC,GERM.N1.wholegene[[3]]$control_CNV)

#####Test messing around with CNV count corrections
#Read test data
test.dat <- SEIZ.N1.wholegene[[3]]
#Scale log-transformed gene length, log-transformed exonic bases, and GC content
test.dat$GC.norm <- scale(test.dat$GC,center=T,scale=T)
test.dat$gene_length.norm <- scale(log(test.dat$gene_length),center=T,scale=T)
test.dat$exonic_bases.norm <- scale(log(test.dat$exonic_bases),center=T,scale=T)
#Fit of normalized length, GC, and exonic bases - all CNVs
triple_norm.fit <- lm(all_CNV ~ GC.norm + gene_length.norm + exonic_bases.norm,data=test.dat)
summary(triple_norm.fit)
#Correct for normalized GC content, length, and exonic bases - all CNVs
adjustments <- ((triple_norm.fit$coefficients[2]*test.dat$GC.norm)+
                  (triple_norm.fit$coefficients[3]*test.dat$gene_length.norm)+
                  (triple_norm.fit$coefficients[4]*test.dat$exonic_bases.norm))
test.dat$all_CNV.trip_norm_adj <- round(test.dat$all_CNV-adjustments)
test.dat$all_CNV.trip_norm_adj[which(test.dat$all_CNV.trip_norm_adj<0)] <- 0
hist(test.dat$all_CNV,breaks=100,col="black",ylim=c(0,9000))
hist(test.dat$all_CNV.trip_norm_adj,breaks=100,col="black",ylim=c(0,9000))
plot(test.dat$all_CNV,test.dat$all_CNV.trip_norm_adj)
#Triple correction for cases alone
triple_norm.cases.fit <- lm(case_CNV ~ GC.norm + gene_length.norm + exonic_bases.norm,data=test.dat)
summary(triple_norm.cases.fit)
adjustments <- ((triple_norm.cases.fit$coefficients[2]*test.dat$GC.norm)+
                  (triple_norm.cases.fit$coefficients[3]*test.dat$gene_length.norm)+
                  (triple_norm.cases.fit$coefficients[4]*test.dat$exonic_bases.norm))
test.dat$case_CNV.trip_norm_adj <- round(test.dat$case_CNV-adjustments)
test.dat$case_CNV.trip_norm_adj[which(test.dat$case_CNV.trip_norm_adj<0)] <- 0
hist(test.dat$case_CNV,breaks=100,col="black",ylim=c(0,9000))
hist(test.dat$case_CNV.trip_norm_adj,breaks=100,col="black",ylim=c(0,9000))
plot(test.dat$case_CNV,test.dat$case_CNV.trip_norm_adj)
#Triple correction for controls alone
triple_norm.controls.fit <- lm(control_CNV ~ GC.norm + gene_length.norm + exonic_bases.norm,data=test.dat)
summary(triple_norm.controls.fit)
adjustments <- ((triple_norm.controls.fit$coefficients[2]*test.dat$GC.norm)+
                  (triple_norm.controls.fit$coefficients[3]*test.dat$gene_length.norm)+
                  (triple_norm.controls.fit$coefficients[4]*test.dat$exonic_bases.norm))
test.dat$control_CNV.trip_norm_adj <- round(test.dat$control_CNV-adjustments)
test.dat$control_CNV.trip_norm_adj[which(test.dat$control_CNV.trip_norm_adj<0)] <- 0
hist(test.dat$control_CNV,breaks=100,col="black",ylim=c(0,9000))
hist(test.dat$control_CNV.trip_norm_adj,breaks=100,col="black",ylim=c(0,9000))
plot(test.dat$control_CNV,test.dat$control_CNV.trip_norm_adj)



plot(test.dat$control_CNV.trip_norm_adj/nsamp[1],
     test.dat$case_CNV.trip_norm_adj/nsamp[2])
abline(0,1)
theta <- test.dat$case_CNV.trip_norm_adj/nsamp[2]-test.dat$control_CNV.trip_norm_adj/nsamp[1]
theta.z <- scale(theta,scale=T,center=F)
hist(theta.z,breaks=100,col="black",xlim=c(-8,8))
theta.z.p <- pnorm(theta.z,lower.tail=F)
hist(theta.z.p,breaks=50)
cleanQQ(pvector=theta.z.p)
theta.z.q <- p.adjust(theta.z.p,method="fdr")
#Nominally signif
length(which(theta.z.p<0.05))
test.dat[(which(theta.z.p<0.05)),1]
#FDR-corrected signif
length(which(theta.z.q<0.05))
test.dat[(which(theta.z.q<0.05)),1]
#Bonferroni-corrected signif
length(which(theta.z.p<0.05/length(theta.z.p)))
test.dat[(which(theta.z.p<0.05/length(theta.z.p))),1]

GC.fit <- lm(all_CNV ~ GC.norm,data=test.dat)
plot(test.dat$GC.norm,test.dat$all_CNV,
     pch=19,col=adjustcolor("black",alpha=0.1))
abline(GC.fit,col="red")
#Correct for normalized GC content alone
test.dat$all_CNV.GC_adj <- round(test.dat$all_CNV-((GC.fit$coefficients[2]*test.dat$GC.norm)))
test.dat$all_CNV.GC_adj[which(test.dat$all_CNV.GC_adj<0)] <- 0
hist(test.dat$all_CNV.GC_adj,breaks=100,col="black")

#Gene length alone
length.fit <- lm(all_CNV ~ log(gene_length),data=test.dat)
plot(log(test.dat$gene_length),test.dat$all_CNV)
abline(length.fit,col="red")
test.dat$all_CNV.len_adj <- test.dat$all_CNV-(length.fit$coefficients[2]*log(test.dat$gene_length)+length.fit$coefficients[1])
length.refit <- lm(all_CNV.len_adj ~ log(gene_length),data=test.dat)
plot(log(test.dat$gene_length),test.dat$all_CNV.len_adj)
abline(length.refit,col="red")

#GC alone
GC.fit <- lm(all_CNV ~ GC,data=test.dat)
plot(test.dat$GC,test.dat$all_CNV,pch=19,col=adjustcolor("black",alpha=0.1))
abline(GC.fit,col="red")
test.dat$all_CNV.GC_adj <- round(test.dat$all_CNV-((GC.fit$coefficients[2]*test.dat$GC)+GC.fit$coefficients[1]))
test.dat$all_CNV.GC_adj[which(test.dat$all_CNV.GC_adj<0)] <- 0
GC.refit <- lm(all_CNV.GC_adj ~ GC,data=test.dat)
plot(test.dat$GC,test.dat$all_CNV.GC_adj,pch=19,col=adjustcolor("black",alpha=0.1))
abline(GC.refit,col="red")
plot(test.dat$all_CNV,round(test.dat$all_CNV.GC_adj))
#Length after GC correction
length <- lm(all_CNV.GC_adj ~ log(gene_length),data=test.dat)
plot(log(test.dat$gene_length),test.dat$all_CNV.GC_adj)


#Exonic bases alone
exbases.fit <- lm(all_CNV ~ log(exonic_bases),data=test.dat)
plot(log(test.dat$exonic_bases),test.dat$all_CNV)
abline(exbases.fit,col="red")
#Gene length and exonic bases combined
len_ex.fit <- lm(all_CNV ~ log(gene_length) + log(exonic_bases),data=test.dat)
test.dat$all_CNV.dual_adj <- (test.dat$all_CNV+
                                len_ex.fit$coefficients[1]+
                                (len_ex.fit$coefficients[2]*log(test.dat$gene_length))+
                                   (len_ex.fit$coefficients[3]*log(test.dat$exonic_bases)))
hist(test.dat$all_CNV.dual_adj,breaks=100)
hist(test.dat$all_CNV,breaks=100)
plot(test.dat$all_CNV,round(test.dat$all_CNV.dual_adj))
abline(lm(all_CNV.dual_adj ~ all_CNV, dat=test.dat),col="red")

#####Function to make scaled dotplot
plotDots <- function(df,CNV,nCTRL,nCASE,sim=10000,colP=T){
  #Gather plotting data
  dat <- data.frame(GERM.N1.exonic[,CNV]$control_CNV/nCTRL,
                    GERM.N1.exonic[,CNV]$case_CNV/nCASE)
  names(dat) <- c("CTRL","CASE")
  #Calculate Beta-Binomial ratio tests
  if(colP==T){
    # #CASE > CTRL
    # CASE.theta <- rbeta(sim,
    #                     (dat$CASE)+1,
    #                     (nCASE-dat$CASE)+1)
    # CTRL.theta <- rbeta(sim,
    #                     (dat$CTRL)+1,
    #                     (nCTRL-dat$CTRL)+1)
    # d.theta <- CASE.theta-CTRL.theta
    # thresholds <- quantile(d.theta,c(0.05/nrow(dat),0.05,0.95,1-(0.05/nrow(dat))))
    # #Test Beta-Binomial plot
    # plot(density(d.theta),
    #      xlab="CASE.theta - CTRL.theta",
    #      ylab="p(CASE.theta - CTRL.theta | nCASE, nCTRL)",
    #      main="Bayesian Posterior Probability of case CNV =/= control CNV",
    #      frame.plot=FALSE,cex.lab=1.5,lwd=3,yaxt="no")
    # abline(v=thresholds,lty=c(2,1,1,2),col=c("black","black","forestgreen","forestgreen"))
    # length(which(as.numeric(d.theta)>=as.numeric(thresholds[4])))
    # length(which(as.numeric(d.theta)>=as.numeric(thresholds[3])))
    # length(which(as.numeric(d.theta)<=as.numeric(thresholds[2])))
    # length(which(as.numeric(d.theta)<=as.numeric(thresholds[1])))
    #Fisher's exact apply
    pvals <- sapply(1:500,function(i){
      fisher.test(matrix(c(nCTRL-dat[i,1]*nCTRL,
                           nCASE-dat[i,2]*nCASE,
                           nCTRL,nCASE),
                         byrow=T,nrow=2))$p.value
    })
  }
  #Instantiate color vector
  colvect.sub <- sapply(pvals,function(p){
    if(p<=0.05){
      return("red")
    }else{
      return("black")
    }
  })
  colvect <- c(colvect.sub,rep(NA,nrow(dat)-length(sub)))
  #Set plotting parameters
  par(mar=c(4,4,0.5,0.5))
  #Determine axis maximum
  axMax <- max(dat)
  #Plot
  plot(dat$CTRL,dat$CASE,
       xlim=c(0,axMax),ylim=c(0,axMax),
       xlab="Control CNVs per Sample",
       ylab="Case CNVs per Sample",
       bg=colvect,pch=21)
  # #Overlay density plot
  # complete.dat <- dat[complete.cases(dat),]
  # contours <- kde2d(complete.dat[,1],complete.dat[,2],n=1000)
  # contour(contours,drawlabels=FALSE,nlevels=25,add=TRUE,
  #         col=adjustcolor(rev(rainbow(45)[1:32]),alpha=0.5))
  #X=Y line
  abline(0,1)
  #Linear regression line
  # abline(lm(nCASE ~ nCTRL,data=dat))
}

#####Calculate nominal and Bonferroni 95% CI bounds
# #Nominal
# GERM.nomP <- sapply(0:nsamp[1],function(nCTRL){
#   sapply(0:nsamp[2],function(nGERM){
#     fisher.test(matrix(c(nsamp[1]-nCTRL,nsamp[2]-nCASE,
#                          nCTRL,nCASE),
#                        byrow=T,nrow=2))$p.value
#   })
# })
# #Bonferroni
# GERM.bonP <- GERM.nomP*nrow(GERM[,1])

#####Test plots - case vs control
# #GERM E2 CNV
# plot(GERM.E2.exonic[,1]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,1]$case_CNV/nsamp[2])
# abline(0,1)
# #GERM E2 DEL
# plot(GERM.E2.exonic[,2]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,2]$case_CNV/nsamp[2])
# abline(0,1)
# #GERM E2 DUP
# plot(GERM.E2.exonic[,3]$control_CNV/nsamp[1],
#      GERM.E2.exonic[,3]$case_CNV/nsamp[2])
# abline(0,1)
#GERM N1 exonic CNV
plotDots(GERM.N1.exonic,1,nsamp[1],nsamp[2],colP=F)
#GERM N1 exonic DEL
plotDots(GERM.N1.exonic,2,nsamp[1],nsamp[2],colP=F)
#GERM N1 exonic DUP
plotDots(GERM.N1.exonic,3,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene CNV
plotDots(GERM.N1.wholegene,1,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene DEL
plotDots(GERM.N1.wholegene,2,nsamp[1],nsamp[2],colP=F)
#GERM N1 whole-gene DUP
plotDots(GERM.N1.wholegene,3,nsamp[1],nsamp[2],colP=F)

#####Fisher test
dat <- head(as.data.frame(GERM.N1.exonic[,2]),500)
i <- which(dat$case_CNV==max(dat$case_CNV))
fisher.test(matrix(c(as.numeric(dat$control_CNV[i]),
                     as.numeric(dat$case_CNV[i]),
                     nsamp[1]-as.numeric(dat$control_CNV[i]),
                     nsamp[2]-as.numeric(dat$case_CNV[i])),
                   byrow=T,nrow=2))$p.value

