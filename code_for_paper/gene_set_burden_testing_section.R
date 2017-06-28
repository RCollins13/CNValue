#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to perform analyses required for gene set enrichment section of rCNV paper


###################
#####Set parameters
###################
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
phenos <- c("GERM","UNK","NEURO","NDD","DD","PSYCH","SCZ","ASD","SEIZ",
            "HYPO","BEHAV","ID","SOMA","HEAD","GRO","CARD","SKEL","DRU",
            "MSC","EE","INT","EMI","CNCR","CGEN","CSKN","CGST","CRNL",
            "CBRN","CLNG","CBST","CEND","CHNK","CLIV","CMSK","CBLD")



######################
#####Set color vectors
######################
cols.CTRL <- c("#A5A6A7","#DCDDDF","#EAEBEC","#F8F8F9")
cols.GERM <- c("#7B2AB3","B07FD1","#CAAAE1","#E5D4F0")
cols.NEURO <- c("#00BFF4","#66D9F8","#99E5FB","#CCF2FD")
cols.SOMA <- c("#EC008D","#F466BB","#F799D1","#FBCCE8")
cols.CNCR <- c("#FFCB00","#FFCB00","#FFE066","#FFF5CC")



####################################################
#####Enrichment of protein-coding genes in all genes
####################################################
#Load data
OR <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                            "CNV_E2_all.effectSizes.txt",sep=""),header=F)
names(OR) <- c("anno",phenos)
pvals <- read.table(paste(WRKDIR,"plot_data/allSet_burden_results/",
                           "CNV_E2_all.pvals.txt",sep=""),header=F)
names(pvals) <- c("anno",phenos)
#Scrape ORs & p-vals
OR$GERM[which(OR$anno=="GENESET.All_protein_coding_genes")]
pvals$GERM[which(OR$anno=="GENESET.All_protein_coding_genes")]
OR$CNCR[which(OR$anno=="GENESET.All_protein_coding_genes")]
pvals$CNCR[which(OR$anno=="GENESET.All_protein_coding_genes")]




#######################################
#####Get background rates of enrichment
#######################################
#All gene sets, all phenotypes (exonic)
OR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                       "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
names(OR) <- c("anno",phenos)
OR.means <- apply(OR[,-1],1,mean,na.rm=T)
mean(OR.means,na.rm=T)
options(scipen=-100)
wilcox.test(OR.means,alternative="greater")$p.value
options(scipen=100)
#All gene sets, all phenotypes (whole-gene)
OR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                       "CNV_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
names(OR) <- c("anno",phenos)
mean(apply(OR[,-1],1,mean,na.rm=T),na.rm=T)
#All gene set enrichments, del vs dup, Mann-Whitney (exonic)
DEL <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DEL.means <- apply(DEL[,-1],1,mean,na.rm=T)
DEL.means <- log2(DEL.means)
DUP <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DUP_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DUP.means <- apply(DUP[,-1],1,mean,na.rm=T)
DUP.means <- log2(DUP.means)
wilcox.test(DUP.means,DEL.means)
boxplot(DEL.means,DUP.means,col=c("red","blue"),notch=T)
#All gene set enrichments, del vs dup, Mann-Whitney (whole-gene)
DEL <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DEL.means <- apply(DEL[,-1],1,mean,na.rm=T)
DEL.means <- log2(DEL.means)
DUP <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DUP_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DUP.means <- apply(DUP[,-1],1,mean,na.rm=T)
DUP.means <- log2(DUP.means)
wilcox.test(DUP.means,DEL.means,alternative="greater")
boxplot(DEL.means,DUP.means,col=c("red","blue"),notch=T)
#All gene set enrichments, whole-gene del vs exonic del, Mann-Whitney
DEL.e <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DEL.emeans <- apply(DEL.e[,-1],1,mean,na.rm=T)
DEL.emeans <- log2(DEL.emeans)
DEL.wg <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DEL.wgmeans <- apply(DEL.wg[,-1],1,mean,na.rm=T)
DEL.wgmeans <- log2(DEL.wgmeans)
wilcox.test(DEL.wgmeans,DEL.emeans,alternative="greater")
boxplot(DEL.emeans,DEL.wgmeans,col=c("red","blue"),notch=T)



#################################
#####GERM vs CNCR CNV correlation
#################################
GERM <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
GERM.means <- apply(GERM[,2:23],1,mean,na.rm=T)
GERM.means <- log2(GERM.means)
CNCR <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "CNV_E2_all_exonic.effectSizes.txt",sep=""),header=F)
CNCR.means <- apply(CNCR[,-c(1:23)],1,mean,na.rm=T)
mean(CNCR.means,na.rm=T)
CNCR.means <- log2(CNCR.means)

plot(GERM.means,
     CNCR.means,
     xlim=log2(c(1/5,5)),ylim=log2(c(1/5,5)),
     pch=19,col=adjustcolor("black",alpha=0.3))
abline(h=0,v=0)
abline(lm(CNCR.means ~ GERM.means),lty=2)
cor(GERM.means,
    CNCR.means,
    use="complete.obs",
    method="spearman")
options(scipen=-1000)
cor.test(GERM.means,
         CNCR.means,
         method="spearman")$p.value
options(scipen=1000)



##############################
#####DEL vs DUP OR correlation
##############################
#Exonic
DEL <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "DEL_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DEL.means <- apply(DEL[,-1],1,mean,na.rm=T)
DEL.means <- log2(DEL.means)
DUP <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                         "DUP_E2_all_exonic.effectSizes.txt",sep=""),header=F)
DUP.means <- apply(DUP[,-1],1,mean,na.rm=T)
mean(DUP.means,na.rm=T)
DUP.means <- log2(DUP.means)
plot(DEL.means,
     DUP.means,
     xlim=log2(c(1/5,5)),ylim=log2(c(1/5,5)),
     pch=19,col=adjustcolor("black",alpha=0.3))
abline(h=0,v=0)
abline(lm(DUP.means ~ DEL.means),lty=2)
cor(DEL.means,
    DUP.means,
    use="complete.obs",
    method="spearman")
options(scipen=-1000)
cor.test(DEL.means,
         DUP.means,
         method="spearman")$p.value
options(scipen=1000)
#Exonic
DEL <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DEL_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DEL.means <- apply(DEL[,-1],1,mean,na.rm=T)
DEL.means <- log2(DEL.means)
DUP <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",
                        "DUP_E2_all_wholegene.effectSizes.txt",sep=""),header=F)
DUP.means <- apply(DUP[,-1],1,mean,na.rm=T)
mean(DUP.means,na.rm=T)
DUP.means <- log2(DUP.means)
plot(DEL.means,
     DUP.means,
     xlim=log2(c(1/5,5)),ylim=log2(c(1/5,5)),
     pch=19,col=adjustcolor("black",alpha=0.3))
abline(h=0,v=0)
abline(lm(DUP.means ~ DEL.means),lty=2)
cor(DEL.means,
    DUP.means,
    use="complete.obs",
    method="spearman")
options(scipen=-1000)
cor.test(DEL.means,
         DUP.means,
         method="spearman")$p.value
options(scipen=1000)



#######################################
#####SIGNIFICANT SET BREAKDOWN BY PHENO
#######################################
#Read p-values
pvals <- lapply(c("CNV","DEL","DUP"),function(CNV){
  dat <- read.table(paste(WRKDIR,"plot_data/geneSet_burden_results/",CNV,
                          "_E2_all_exonic.pvals.txt",sep=""),header=F)
  dat[,-1] <- apply(dat[,-1],2,as.numeric)
  colnames(dat) <- c("geneSet",phenos)
  return(dat)
})
#Iterate over CNV types and count # signif phenotypes
signif <- lapply(pvals,function(dat){
  apply(dat[,-1],1,function(vals){
    return(length(which(vals<=0.05/263 & !is.na(vals))))
  })
})
#Count significant sets by 0, 1, 2-34, 35
signif.counts <- lapply(signif,function(vals){
  none <- length(which(vals==0))
  notNone <- length(which(vals>0))
  one <- length(which(vals==1))
  some <- length(which(vals>1 & vals<35))
  all <- length(which(vals==35))
  total <- sum(c(none,one,some,all))
  return(data.frame(c(none,notNone,one,some,all,total),
                    c(none,notNone,one,some,all,total)/length(vals)))
})
#Count number of significant sets per pheno group
pheno.signif <- lapply(pvals,function(dat){
  apply(dat[,-1],2,function(vals){
    return(length(which(vals<=0.05/263 & !is.na(vals))))
  })
})
#Median per cancer group and per germline group
median(pheno.signif[[1]][1:22])
median(pheno.signif[[1]][-c(1:22)])
options(scipen=-100)
wilcox.test(pheno.signif[[1]][-c(1:22)],
            pheno.signif[[1]][1:22])$p.value
#Del vs dup for germline
wilcox.test(pheno.signif[[2]][1:22],
            pheno.signif[[3]][1:22])
#CNCR vs all other cancer types
wilcox.test(pheno.signif[[1]][-c(1:23)],
            mu=pheno.signif[[1]][23])
#GERM vs all other germline groups
wilcox.test(pheno.signif[[1]][2:22],
            mu=pheno.signif[[1]][1])
options(scipen=1000)

#Test plot
plot(c(1,35),c(0,100),type="n")
points(pheno.signif[[1]],pch=19,col="gray20")
points(pheno.signif[[2]],pch=19,col="red")
points(pheno.signif[[3]],pch=19,col="blue")


############################################################
#####Analysis of private and universal significant gene sets
############################################################
#Match gene sets private to single phenotypes with their phenos - CNV
data.frame(which(signif[[1]]==1),
           DEL[which(signif[[1]]==1),1],
           sapply(which(signif[[1]]==1),function(i){
             phenos[which(pvals[[1]][i,-1]<=0.05/nrow(pvals[[1]]))]
           }))
#LCL genes in blood cancer - CNV
plot(-log10(as.numeric(pvals[[1]][78,-1])))
#cartilage development in integument defects - CNV
plot(-log10(as.numeric(pvals[[1]][240,-1])))
#head/neck in ID - CNV
plot(-log10(as.numeric(pvals[[1]][108,-1])))
#Match gene sets private to single phenotypes with their phenos - DEL
data.frame(which(signif[[2]]==1),
           DEL[which(signif[[2]]==1),1],
           sapply(which(signif[[2]]==1),function(i){
             phenos[which(pvals[[2]][i,-1]<=0.05/nrow(pvals[[2]]))]
           }))
#ASD genes in ASD - DEL
plot(-log10(as.numeric(pvals[[2]][21,-1])))
#Hematopoietic cell diff genes in liver cancer - DEL
plot(-log10(as.numeric(pvals[[2]][239,-1])))
#Match gene sets private to single phenotypes with their phenos - DUP
data.frame(which(signif[[3]]==1),
           DEL[which(signif[[3]]==1),1],
           sapply(which(signif[[3]]==1),function(i){
             phenos[which(pvals[[3]][i,-1]<=0.05/nrow(pvals[[3]]))]
           }))
#Nervous system dev genes in seizures - DUP
plot(-log10(as.numeric(pvals[[3]][228,-1])))
#Write list of gene sets significant to all, one, and no phenos
write.table(DEL[which(signif[[1]]==35),1],
            paste(WRKDIR,"/rCNVmap/misc/geneSets_signif_allPhenos.list",sep=""),
            col.names=F,row.names=F,quote=F)
write.table(DEL[which(signif[[1]]==1),1],
            paste(WRKDIR,"/rCNVmap/misc/geneSets_signif_onePheno.list",sep=""),
            col.names=F,row.names=F,quote=F)
write.table(DEL[which(signif[[1]]==0),1],
            paste(WRKDIR,"/rCNVmap/misc/geneSets_signif_noPhenos.list",sep=""),
            col.names=F,row.names=F,quote=F)
#Match gene sets unique to N phenotypes with their phenos - CNV
N=5 #number of phenotypes
i=3 #CNV type
data.frame(which(signif[[i]]==N),
           DEL[which(signif[[i]]==N),1],
           sapply(which(signif[[i]]==N),function(j){
             paste(unlist(phenos[which(pvals[[i]][j,-1]<=0.05/nrow(pvals[[i]]))]),collapse=",")
           }))










