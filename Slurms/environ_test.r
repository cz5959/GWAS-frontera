#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
setwd("/scratch1/08005/cz5959/QC/MFI/autosome")

# parameters
option_list = list(
  make_option(c("-i","--snps"), type="integer", default=100,help="i: num of causal SNPs",metavar="character"), 
  make_option(c("-g","--heritability"), type="numeric", default=100,help="g: h2",metavar="character"), 
  make_option(c("-e","--environment"), type="numeric", default=100,help="e: environment",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

snp_num <- opt$snps; print(snp_num)
h2 <- opt$heritability; print(h2)
E_ratio <- opt$environment; print(E_ratio)

# create sample of individuals with random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

# set seed
set.seed(1)

# genotype matrix for 300k individuals
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; snp_freqs <- snp_freqs$x    # load 20k total snp frequencies
setwd("/scratch1/08005/cz5959/simulation")
load("simulation_matrix_k_5k.RData")  # load matrix of 300k by 5k snp genotypes 

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]   # subset i causal snp frequencies
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]  # subset i causal snp genotypes

######################################################################################

# Beta (effect size)
Beta <- rnorm(snp_num, mean=0, sd=1)    # get vector of i Beta effect sizes

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G <- (Beta^2) * snp_freqs * (1-snp_freqs)   # genetic variance
  var_E_m <- (var_G * (1-h2)) / h2          # environmental variance (males)
  var_E_f <- var_E_m * E_ratio              # environmental variance for female based on proportion
  E <- NULL
  for (i in 1:length(sex)) {    # assign environmental effect based on gender
    if (sex[i] == 0) {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_f))
    } else {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_m))
    }
  }
  return(E)
}

# environment and heritability combination
E <- get_environment(h2, E_ratio)   # call function to get vector of environmental effects
rm(snp_freqs)
pheno <- rowSums(t(t(genotype_matrix_i)*Beta)) + E    # get vector of 300k phenotypes
rm(genotype_matrix_i)

# GWAS function for each SNP
GWAS <- function(genotype, pheno, E) {
  # linear regression of phenotype on genotype + environmental effect
  f_model <- lm(pheno[f_id] ~ genotype[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ genotype[m_id] + E[m_id])
  # grab Beta and SE coefficients for the SNP
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

# k snps, LD blocks
mash_BETA <- NULL
mash_SE <- NULL
for (j in c(1,10,15,20)) {    # for each of the premade genotype matrices
  if (j != 1){
    load(paste0("simulation_matrix_k_",j,"k.RData"))    # go through each SNP
  }
  for (i in 1:length(colnames(genotype_matrix_k))) {
    gwas <- GWAS(genotype_matrix_k[,i], pheno, E)     # call GWAS function
    mash_BETA <- rbind(mash_BETA, data.frame(female= gwas[1,1], male= gwas[2,1]))
    mash_SE <- rbind(mash_SE, data.frame(female= gwas[1,2], male= gwas[2,2]))
  }
  print(j)
  save(mash_BETA, mash_SE, file="simulation_mash_wip.RData")
  rm(genotype_matrix_k)
}

save(mash_BETA, mash_SE, file=paste0("mash_",snp_num,"_",h2,"_",E_ratio,".RData"))

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
write.table(mixture, file=paste0(snp_num,"_",h2,"_",E_ratio,".txt"),row.names = FALSE, sep="\t")
