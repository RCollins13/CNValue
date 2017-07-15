#!/usr/bin/env Rscript

##############################
#    Rare CNV Map Project    #
##############################

# Copyright (c) 2017 Ryan L. Collins
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits availble on GitHub

#Script to run per-gene CNV burden model

##Command line arguments:
# infile: path to geneScore input file
# nCTRL: number of control samples
# nCASE: number of case samples
##Command line option:
# outfile: path to output file

#####Set params
options(scipen=1000,stringsAsFactors=F)

#####Load requirements
require("optparse")

##Local dev testing parameters
# infile <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/plot_data/perGene_burden/NDD_DEL_E4_exonic.geneScore_data.txt"
# #infile <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/plot_data/perGene_burden/CNCR_DEL_E2_exonic.geneScore_data.txt"
# #infile <- "/Users/rlc/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/plot_data/perGene_burden/NDD_DUP_E4_wholegene.geneScore_data.txt"
# infile <- "~/scratch/ASD_DEL_E4_exonic.geneScore_data.txt"
# outfile <- "~/scratch/geneScore.test.output.txt"
# nCTRL <- 38628
# nCASE <- 3398
#
# #Test run
# dat <- readGeneScores(infile)
# df <- adjustCounts(dat)
# fisher.results <- calcFisherStats(df,nCTRL,nCASE)
# ratio.results <- calcRatioStats(df,nCTRL,nCASE)
# all.results <- cbind(fisher.results,ratio.results[,-c(1:25)])

######################################################
#####Helper function to read & clean geneScore dataset
######################################################
readGeneScores <- function(path){
  #Read file
  df <- read.table(path,header=T,comment.char="",sep="\t")

  #Clean file
  df <- as.data.frame(df)
  names(df)[1] <- "gene"
  df$all_CNV <- df$case_CNV+df$control_CNV
  df$all_CNV_weighted <- df$case_CNV_weighted+df$control_CNV_weighted

  #Return file
  return(df)
}

########################################################
#####Helper function to adjust CNV counts by GC content,
#####gene length, and number of exonic bases
########################################################
adjustCounts <- function(df){
  #Normalize log-transformed gene length, log-transformed exonic bases, and GC content
  df$GC.norm <- scale(df$GC,center=T,scale=T)
  df$gene_length.norm <- scale(log(df$gene_length),center=T,scale=T)
  df$exonic_bases.norm <- scale(log(df$exonic_bases),center=T,scale=T)

  #Linear fit of raw CNV counts by normalized length, GC, and exonic bases
  triple_norm.CTRL.fit <- lm(control_CNV ~ GC.norm + gene_length.norm + exonic_bases.norm,data=df)
  triple_norm.CASE.fit <- lm(case_CNV ~ GC.norm + gene_length.norm + exonic_bases.norm,data=df)

  #Calculate corrections for raw CNV counts by normalized GC content, length, and exonic bases
  adjustments.CTRL <- ((triple_norm.CTRL.fit$coefficients[2]*df$GC.norm)+
                         (triple_norm.CTRL.fit$coefficients[3]*df$gene_length.norm)+
                         (triple_norm.CTRL.fit$coefficients[4]*df$exonic_bases.norm))
  adjustments.CASE <- ((triple_norm.CASE.fit$coefficients[2]*df$GC.norm)+
                         (triple_norm.CASE.fit$coefficients[3]*df$gene_length.norm)+
                         (triple_norm.CASE.fit$coefficients[4]*df$exonic_bases.norm))

  #Apply corrections for raw CNV counts and round negative corrected counts to zero
  df$control_CNV.adj <- df$control_CNV-adjustments.CTRL
  df$control_CNV.adj_round <- round(df$control_CNV.adj)
  df$control_CNV.adj_round[which(df$control_CNV.adj_round<0)] <- 0
  df$case_CNV.adj <- df$case_CNV-adjustments.CASE
  df$case_CNV.adj_round <- round(df$case_CNV.adj)
  df$case_CNV.adj_round[which(df$case_CNV.adj_round<0)] <- 0

  #Sum adjusted case+control CNV counts
  df$sum_adj_CNV <- df$control_CNV.adj + df$case_CNV.adj
  df$sum_adj_CNV_round <- df$control_CNV.adj_round + df$case_CNV.adj_round

  #Linear fit of weighted CNV counts by normalized length, GC, and exonic bases
  triple_norm.CTRL_weighted.fit <- lm(control_CNV_weighted ~ GC.norm + gene_length.norm + exonic_bases.norm,data=df)
  triple_norm.CASE_weighted.fit <- lm(case_CNV_weighted ~ GC.norm + gene_length.norm + exonic_bases.norm,data=df)

  #Calculate corrections for raw CNV counts by normalized GC content, length, and exonic bases
  adjustments.CTRL_weighted <- ((triple_norm.CTRL_weighted.fit$coefficients[2]*df$GC.norm)+
                         (triple_norm.CTRL_weighted.fit$coefficients[3]*df$gene_length.norm)+
                         (triple_norm.CTRL_weighted.fit$coefficients[4]*df$exonic_bases.norm))
  adjustments.CASE_weighted <- ((triple_norm.CASE_weighted.fit$coefficients[2]*df$GC.norm)+
                         (triple_norm.CASE_weighted.fit$coefficients[3]*df$gene_length.norm)+
                         (triple_norm.CASE_weighted.fit$coefficients[4]*df$exonic_bases.norm))

  #Apply corrections for raw CNV counts and round negative corrected counts to zero
  df$control_CNV_weighted.adj <- df$control_CNV_weighted-adjustments.CTRL_weighted
  df$control_CNV_weighted.adj_round <- round(df$control_CNV_weighted.adj)
  df$control_CNV_weighted.adj_round[which(df$control_CNV_weighted.adj_round<0)] <- 0
  df$case_CNV_weighted.adj <- df$case_CNV_weighted-adjustments.CASE_weighted
  df$case_CNV_weighted.adj_round <- round(df$case_CNV_weighted.adj)
  df$case_CNV_weighted.adj_round[which(df$case_CNV_weighted.adj_round<0)] <- 0

  #Sum adjusted case+control CNV counts
  df$sum_adj_CNV_weighted <- df$control_CNV_weighted.adj + df$case_CNV_weighted.adj
  df$sum_adj_CNV_weighted_round <- df$control_CNV_weighted.adj_round + df$case_CNV_weighted.adj_round

  #Reorder & reheader columns
  df.out <- df[,c(1,2,12,3,13,4,11,
                  6,14,15,8,20,21,
                  5,16,17,7,22,23,
                  9,18,19,10,24,25)]
  names(df.out) <- c("gene","gene_length","gene_length_Z","exonic_bases","exonic_bases_Z","GC","GC_Z",
                     "control_raw_CNV","control_raw_CNV_adjusted","control_raw_CNV_adjusted_rounded",
                     "control_weighted_CNV","control_weighted_CNV_adjusted","control_weighted_CNV_adjusted_rounded",
                     "case_raw_CNV","case_raw_CNV_adjusted","case_raw_CNV_adjusted_rounded",
                     "case_weighted_CNV","case_weighted_CNV_adjusted","case_weighted_CNV_adjusted_rounded",
                     "total_raw_CNV","total_raw_CNV_adjusted","total_raw_CNV_adjusted_rounded",
                     "total_weighted_CNV","total_weighted_CNV_adjusted","total_weighted_CNV_adjusted_rounded")

  #Return data frame
  return(df.out)
}

###############################################################
#####Helper function to calculate Fisher's Exact ORs and p-vals
###############################################################
#Note: uses rounded raw CNV counts, not weighted counts
calcFisherStats <- function(df,nCTRL,nCASE){
  #Get unique pairs of case/control CNV counts
  countsTable <- table(df$control_raw_CNV_adjusted_rounded,
                       df$case_raw_CNV_adjusted_rounded)
  pairs <- sapply(1:nrow(countsTable),function(r){
    rowval <- as.numeric(rownames(countsTable)[r])
    matches <- as.numeric(colnames(countsTable)[which(countsTable[r,]>0)])
    pairs <- t(data.frame(rep(rowval,times=length(matches)),matches))
    return(pairs)
  })
  pairs <- as.data.frame(matrix(unlist(pairs),ncol=2,byrow=T))
  names(pairs) <- c("control_CNV","case_CNV")

  #Create Fisher Test lookup table for all gnees
  fisher.lookup <- data.frame("control_noCNV"=nCTRL-pairs$control_CNV,
                              "control_CNV"=pairs$control_CNV,
                              "case_noCNV"=nCASE-pairs$case_CNV,
                              "case_CNV"=pairs$case_CNV)
  fisher.results <- t(apply(fisher.lookup,1,function(counts){
    #Case > Control
    case.res <- fisher.test(matrix(counts,nrow=2,byrow=F),alternative="greater")
    #Control > Case
    control.res <- fisher.test(matrix(counts,nrow=2,byrow=F),alternative="less")
    #Return vector of relevant values
    out <- as.numeric(c(case.res$estimate,
                        case.res$conf.int[1],
                        case.res$conf.int[2],
                        case.res$p.value,
                        control.res$estimate,
                        control.res$conf.int[1],
                        control.res$conf.int[2],
                        control.res$p.value))
    return(out)
  }))
  fisher.lookup <- cbind(fisher.lookup,fisher.results)
  names(fisher.lookup)[-c(1:4)] <- c("case_gt_control_Fisher_OR",
                                     "case_gt_control_Fisher_lowerCI",
                                     "case_gt_control_Fisher_upperCI",
                                     "case_gt_control_Fisher_unadjusted_p",
                                     "control_gt_case_Fisher_OR",
                                     "control_gt_case_Fisher_lowerCI",
                                     "control_gt_case_Fisher_upperCI",
                                     "control_gt_case_Fisher_unadjusted_p")

  #Add FDR-corrected q-values and Bonferroni-corrected p-values to Fisher lookup table
  fisher.lookup$case_gt_control_Fisher_FDR_q <- p.adjust(fisher.lookup$case_gt_control_Fisher_unadjusted_p,method="fdr")
  fisher.lookup$case_gt_control_Fisher_Bonferroni_p <- p.adjust(fisher.lookup$case_gt_control_Fisher_unadjusted_p,method="bonferroni")
  fisher.lookup$control_gt_case_Fisher_FDR_q <- p.adjust(fisher.lookup$control_gt_case_Fisher_unadjusted_p,method="fdr")
  fisher.lookup$control_gt_case_Fisher_Bonferroni_p <- p.adjust(fisher.lookup$control_gt_case_Fisher_unadjusted_p,method="bonferroni")

  #Reorder columns for Fisher lookup table
  fisher.lookup <- fisher.lookup[,c(2,4:8,13:14,9:12,15:16)]

  #Iterate over all adjusted rounded CNV counts and pull Fisher stats
  fisher.matched <- as.data.frame(t(apply(data.frame(df$control_raw_CNV_adjusted_rounded,
                   df$case_raw_CNV_adjusted_rounded),
        1,function(vals){
          idx <- which(fisher.lookup$control_CNV==vals[1] &
                         fisher.lookup$case_CNV==vals[2])
          return(as.numeric(fisher.lookup[idx,-c(1,2)]))
        })))
  names(fisher.matched) <- names(fisher.lookup)[-c(1,2)]

  #Merge CNV count data with Fisher results
  results <- cbind(df,fisher.matched)

  #Return results
  return(results)
}

#######################################################################
#####Helper function to calculate Z-scores and p-values from ratio test
#######################################################################
#Note: uses adjusted, non-rounded weighted CNV counts, not raw CNV counts
calcRatioStats <- function(df,nCTRL,nCASE){
  #Calculate CNV ratios
  control_ratio <- df$control_weighted_CNV_adjusted/nCTRL
  case_ratio <- df$case_weighted_CNV_adjusted/nCASE

  #Calculate theta (difference in ratios) & Z-scores
  theta <- case_ratio-control_ratio
  theta.z <- scale(theta,scale=T,center=F)

  #Assign percentile based on theta.z
  theta.z.rank <- rank(theta.z)/length(theta.z)

  #Calculate p-values
  theta.z.p.case_gt_control <- pnorm(theta.z,lower.tail=F)
  theta.z.p.control_gt_case <- pnorm(theta.z,lower.tail=T)

  #Adjust p-values (FDR & Bonferroni)
  theta.z.FDR_q.case_gt_control <- p.adjust(theta.z.p.case_gt_control,method="fdr")
  theta.z.Bonferroni_p.case_gt_control <- p.adjust(theta.z.p.case_gt_control,method="bonferroni")
  theta.z.FDR_q.control_gt_case <- p.adjust(theta.z.p.control_gt_case,method="fdr")
  theta.z.Bonferroni_p.control_gt_case <- p.adjust(theta.z.p.control_gt_case,method="bonferroni")

  #Create data frame of results
  res <- data.frame("control_ratio"=control_ratio,
                    "case_ratio"=case_ratio,
                    "theta"=theta,
                    "theta_Zscore"=theta.z,
                    "theta_centile"=theta.z.rank,
                    "case_gt_control_ratio_uncorrected_p"=theta.z.p.case_gt_control,
                    "case_gt_control_ratio_FDR_q"=theta.z.FDR_q.case_gt_control,
                    "case_gt_control_ratio_Bonferroni_p"=theta.z.Bonferroni_p.case_gt_control,
                    "control_gt_case_ratio_uncorrected_p"=theta.z.p.control_gt_case,
                    "control_gt_case_ratio_FDR_q"=theta.z.FDR_q.control_gt_case,
                    "control_gt_case_ratio_Bonferroni_p"=theta.z.Bonferroni_p.control_gt_case)
  results <- cbind(df,res)

  #Return results
  return(results)
}

#######################################
#####Rscript command line functionality
#######################################
#List of Rscript options
option_list <- list(
  make_option(c("-o", "--outfile"), type="character", default="/dev/stdout",
              help="write output to file [default stdout]",
              metavar="character")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] geneScoreData nCTRL nCASE",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#Assigns positional arguments to named variables
infile <- args$args[1]
nCTRL <- as.numeric(args$args[2])
nCASE <- as.numeric(args$args[3])

#Checks for appropriate positional arguments
if(length(args$args) != 3) {
  print_help(parser)
  stop("Must provide an input geneScore data file and control/case sample sizes")
}

#Reads data
df <- readGeneScores(infile)

#Adjusts CNV counts
df.adj <- adjustCounts(df)

#Runs Fisher tests
fisher.results <- calcFisherStats(df.adj,nCTRL,nCASE)

#Runs ratio tests
ratio.results <- calcRatioStats(df.adj,nCTRL,nCASE)

#Combines Fisher & ratio test statistics
all.results <- cbind(fisher.results,ratio.results[,-c(1:25)])

#Fix header
colnames(all.results)[1] <- "#gene"

#Write out
if(opts$outfile != "/dev/stdout"){
  write.table(all.results,opts$outfile,col.names=T,row.names=F,sep="\t",quote=F)
}else{
  print(all.results,row.names=F)
}
