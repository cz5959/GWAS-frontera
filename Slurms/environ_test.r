#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
setwd("/scratch1/08005/cz5959/QC/MFI/autosome")

# parameters
option_list = list(
  make_option(c("-i","--snps"), type="integer", default=100,help="i: num of causal SNPs",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
snp_num <- opt$snps; print(snp_num)

# create sample of individuals -- random sex
# 0=female  1=male
sex <- rbinom(100000,1,0.5)
f_id <- which(sex==0)
m_id <- which(sex==1)

### SNP ###
# get genotype from binom sampling at site k with allele freq
allele_freq <- read.csv("maf_sample_qc.txt", colClasses = 'numeric')
set.seed(1)

# sample i causal SNPs from k genotypes
snpids <- sample(1:nrow(allele_freq),snp_num)
snp_freqs <- allele_freq$maf[snpids]

# get genotype at site i for all individuals
get_genotypes <- function(snp_freq) {
  allele1 <- rbinom(n=length(sex), size=1, prob=snp_freq)
  allele2 <- rbinom(n=length(sex), size=1, prob=snp_freq)
  genotypes <- allele1 + allele2
  if (length(unique(genotypes)) ==1) {
    return(c(100))
  } else {
    return(genotypes)
  }
}

# get environmental effect
get_environment <- function() {
  var_G <- (Beta^2) * snp_af * (1-snp_af)   # genetic variance
  var_E_m <- (var_G * (1-h2)) / h2          # environmental variance (males)
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

GWAS <- function(E) {
  pheno <- (Beta*snp_count) + E      # calculate phenotypes from Beta, snp i, and E
  f_model <- lm(pheno[f_id] ~ snp_count[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ snp_count[m_id] + E[m_id])
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

# get genotypes and betas for all snps
monomorphs <- NULL
mash_BETA <- NULL
mash_SE <- NULL
for (i in 1:length(allele_freq)) {
  snp_count <- get_genotypes(allele_freq$maf[i])    # get allele count for all individuals\
  snp_af <- allele_freq$maf[i]
  if (length(snp_count) == 1) {       # check if snp is monomorphic
    monomorphs[i] <- snpids[i]
    next
  }
  Beta <- rnorm(1, mean=0, sd=1)          # sample effect size B from normal distribution

  for (h2 in c(0.5, 0.05)) {
    for (E_ratio in c(1, 1.5, 5)) {
      E <- get_environment()
      gwas <- GWAS(E)
      mash_BETA <- rbind(mash_BETA, data.frame(h2= h2, E_ratio=E_ratio, female= gwas[1,1], male= gwas[2,1]))
      mash_SE <- rbind(mash_SE, data.frame(h2= h2, E_ratio=E_ratio, female= gwas[1,2], male= gwas[2,2]))
    }
  }

  # print progress
  if (i %% 10 == 0){
    print(i)
  }
}

print(head(mash_BETA,20))
print(head(mash_SE,20))

mash_func <- function() {
  for (h2 in c(0.5, 0.05)) {
    for (E_ratio in c(1, 1.5, 5)) {
      BETA <- mash_BETA[mash_BETA$h2 == h2 & mash_BETA$E_ratio == E_ratio, c(3,4)]
      SE <- mash_SE[mash_SE$h2 == h2 & mash_SE$E_ratio == E_ratio, c(3,4)]
      BETA <- as.matrix(BETA) ; SE <- as.matrix(SE)
      data <- mash_set_data(BETA, SE, zero_Shat_reset =  2.5e-16)
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
    }
  }
}


setwd("/scratch1/08005/cz5959/simulation")
if (length(monomorphs) == length(snpids)) {
  print("All SNPs are monomorphic :(")
} else {
  # run the mash function
  mash_func()
}

print("done")
