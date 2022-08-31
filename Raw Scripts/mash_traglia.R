#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")

option_list = list(
  make_option(c("-n","--name"), type="character", default="",help="file name",metavar="character"),
  make_option(c("-p","--pheno"), type="character", default="",help="phenotype",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
name <- opt$name; print(name)
pheno <- opt$pheno; print(pheno)

# load in RData 
setwd(paste0("/scratch/08005/cz5959/simulation/1000_",pheno))
load(paste0(name,".RData"))

# traglia
z <- (mash_BETA$female - mash_BETA$male) / sqrt( mash_SE$male^2 + mash_SE$female^2)
p <- 2*pnorm(abs(z), 0, 1, lower.tail=FALSE)
p_id <- which(p<5e-8)

mash_func <- function(mash_BETA, mash_SE) {
  ## mash
  mash_BETA <- as.matrix(mash_BETA) ; mash_SE <- as.matrix(mash_SE)
  data <- mash_set_data(mash_BETA, mash_SE, zero_Shat_reset = .Machine$double.eps)
  # set up covariance matrices
  U.c = cov_canonical(data)
  corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
  effect = c(1.5,2,3)
  for (c in corr) {
  for (e in effect) {
      U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
      U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
  }
  }
  U.c[['equal_-0.25_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
  U.c[['equal_-0.5_1']] <- matrix(c(1,-0.5,-0.5,1),2,2)
  U.c[['equal_-0.75_1']] <- matrix(c(1,-0.75,-0.75,1),2,2)
  U.c[['equal_-1_1']] <- matrix(c(1,-1,-1,1),2,2)
  names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")

  # fit mash model
  m.c <- mash(data,U.c)
  mixture <- get_estimated_pi(m.c)
}
# normal
mixture_normal <- mash_func(mash_BETA, mash_SE)
write.table(mixture_normal, file=paste0(name,".txt"),row.names = FALSE, sep="\t")
# ttest
mash_BETA <- mash_BETA[p_id,]
mash_SE <- mash_SE[p_id,]
mixture_ttest <- mash_func(mash_BETA, mash_SE)
write.table(mixture_ttest, file=paste0(name,"_traglia.txt"),row.names = FALSE, sep="\t")



