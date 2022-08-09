#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
library(MASS, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
setwd("/scratch/08005/cz5959/MFI")

# parameters
option_list = list(
  make_option(c("-i","--snps"), type="integer", default=100,help="i: num of causal SNPs",metavar="character"),
  make_option(c("-e","--environment"), type="numeric", default=100,help="e: environment",metavar="character"),
  make_option(c("-a","--amplification"), type="numeric", default=100,help="a: amplification weight",metavar="character")  
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
snp_num <- opt$snps; print(snp_num)
E_ratio <- opt$environment; print(E_ratio)
a_weight <- (opt$amplification)/100; print(amplification)

# create sample of individuals -- random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

# set seed
set.seed(1)

# genotype matrix for 300k
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; snp_freqs <- snp_freqs$x
setwd("/scratch/08005/cz5959/simulation")
load("simulation_matrix_k_5k.RData")

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]

# Beta (effect size)
sigma_1 <- matrix(c(1,1,1,1), 2, 2) 
sigma_2 <- matrix(c(2,1,1,1), 2, 2) 1
mu <- c(0,0)
Beta <- mvrnorm(snp_num*(1-a_weight), mu=mu, Sigma=sigma_1)
Beta <- rbind(Beta, mvrnorm(snp_num*a_weight, mu=mu, Sigma=sigma_2))

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G_m <- (Beta[,2]^2) * 2 * snp_freqs * (1-snp_freqs)   # genetic variance
  var_E_m <- (var_G_m * (1-h2)) / h2          # environmental variance (males)
  var_E_f <- var_E_m * E_ratio              # environmental variance for female based on proportion
  E <- NULL
  for (i in 1:length(sex)) {
    if (sex[i] == 0) {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_f))
    } else {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_m))
    }
  }
  return(E)
}

# environment and heritability combination
# E_ratio in c(1, 1.5, 5)
E <- get_environment(0.5, E_ratio)
rm(snp_freqs)
genotype_matrix_i <- t(genotype_matrix_i)
genotype_matrix_i[,m_id] <- genotype_matrix_i[,m_id]*Beta[,2]
genotype_matrix_i[,f_id] <- genotype_matrix_i[,f_id]*Beta[,1]
pheno <- colSums(genotype_matrix_i) + E
rm(genotype_matrix_i)


# GWAS ## TODo: more snps per gwas
GWAS <- function(genotype, pheno, E) {
  f_model <- lm(pheno[f_id] ~ genotype[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ genotype[m_id] + E[m_id])
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

# k snps, LD blocks
mash_BETA <- NULL
mash_SE <- NULL
for (j in c(1,10,15,20)) {
  if (j != 1){
    load(paste0("simulation_matrix_k_",j,"k.RData"))
  }
  for (i in 1:length(colnames(genotype_matrix_k))) {
    gwas <- GWAS(genotype_matrix_k[,i], pheno, E)
    mash_BETA <- rbind(mash_BETA, data.frame(female= gwas[1,1], male= gwas[2,1]))
    mash_SE <- rbind(mash_SE, data.frame(female= gwas[1,2], male= gwas[2,2]))
  }
  print(j)
  rm(genotype_matrix_k)
}

save(mash_BETA, mash_SE, file=paste0("mash_",snp_num,"_",a_weight,".RData"))


