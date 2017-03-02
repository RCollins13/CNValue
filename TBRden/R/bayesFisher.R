#!/usr/bin/env R

#Copyright (c) 2016 Ryan Collins and Jake Conway
#Distributed under terms of the MIT License

#bayesFisher: Bayesian alternative to Fisher's Exact test (Beta-Binomial test)

bayesFisher <- function(CTRL.bin,   #control CNVs in test bin
                        CTRL.flank, #control CNVs in flank bin
                        CASE.bin,   #case CNVs in test bin
                        CASE.flank  #case CNVs in flank bin
){
  p <- (CASE.bin+CASE.flank-1)/(CASE.bin+CASE.flank+CTRL.bin+CTRL.flank-2)
}
