#!/usr/bin/env Rscript
libLoc <- "/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/"
library("crayon",lib.loc=libLoc)
library("cli",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)

set.seed(1)
n=300000
get_genotypes <- function(snp_freq) {
  allele1 <- rbinom(n=n, size=1, prob=snp_freq)
  allele2 <- rbinom(n=n, size=1, prob=snp_freq)
  genotypes <- allele1 + allele2
  return(genotypes)
}

setwd("/scratch1/08005/cz5959/QC/MFI/autosome")
gwas_snps <- read.csv("maf_sample_1700e2.txt", colClasses = 'numeric') ; gwas_snps <- gwas_snps$x

setwd("/scratch1/08005/cz5959/simulation")
genotype_matrix_k <- NULL
for (i in 1:length(gwas_snps)) {
  snp <- get_genotypes(gwas_snps[i])    # get allele count for all individuals
  genotype_matrix_k <- suppressMessages(bind_cols(genotype_matrix_k, snp))
  if (i %% 1000 == 0){
    print(i)
    write.table(genotype_matrix_k, file=paste0("simulation_matrix_k_",i,".txt"), sep="\t", row.names=FALSE)
  }
}

setwd("/scratch1/08005/cz5959/simulation")
write.table(genotype_matrix_k, file="simulation_matrix_k.txt", sep="\t", row.names=FALSE)
save(genotype_matrix_k, file="simulation_matrix_k.RData")