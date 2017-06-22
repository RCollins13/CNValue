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
mean(apply(OR[,-1],1,mean,na.rm=T),na.rm=T)
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

