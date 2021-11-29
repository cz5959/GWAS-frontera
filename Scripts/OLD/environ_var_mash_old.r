#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

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

##########################################################################
setwd("/scratch1/08005/cz5959/QC/MFI/autosome")

### SNP ###
# get genotype from binom sampling at site k with allele freq
allele_freq <- read.csv("maf_sample_qc.txt", colClasses = 'numeric')
set.seed(1)

n=100000
genotype_matrix_k <- NULL
genotype_matrix_k <- matrix(nrow=n, ncol=nrow(allele_freq))
for (i in 1:nrow(allele_freq)){
    allele1 <- rbinom(n=n, size=1, prob=allele_freq$maf[i])
    allele2 <- rbinom(n=n, size=1, prob=allele_freq$maf[i])
    genotypes <- allele1 + allele2
    genotype_matrix_k[,i] <- genotypes
    #genotype_matrix_k <- cbind(genotype_matrix_k, genotypes)
    # print progress
    if (i %% 1000 == 0){
        print(i)
    }
}

# remove monomorphism
poly <- which(apply(genotype_matrix_k,2,function(x) length(unique(x))) != 1)
poly_num <- length(poly)
print(paste0("# monomorphisms: ", (nrow(allele_freq)-poly_num)))
genotype_matrix_k <- genotype_matrix_k[,poly]
allele_freq <- allele_freq[poly,]

setwd("/scratch1/08005/cz5959/simulation")
save(genotype_matrix_k, allele_freq, file="simulation_matrix.RData")

##########################################################################
setwd("/scratch1/08005/cz5959/simulation")
load(file="simulation_matrix.RData")

# sample i causal SNPs from k genotypes
snpids <- sample(1:nrow(allele_freq),snp_num)
snp_freqs <- allele_freq$maf[snpids]
snps <- genotype_matrix_k[,snpids]

# sample effect size B from normal distribution
Beta <- rnorm(snp_num, mean=0, sd=1)

### ENV VARIANCE ###
# genetic variance for each individuals
var_G <- colSums((Beta^2) * snp_freqs * (1-snp_freqs))

# environmental variance for males from heritability
get_varE_m <- function(h2) {
    var_E_m <- (var_G * (1-h2)) / h2
    return(var_E_m)
}

# environmental variance for female based on proportion
get_varE_f <- function(E_ratio, var_E_m) {
    var_E_f <- var_E_m * E_ratio
    return(var_E_f)
}

# environmental effect for males and females from normal distribution
get_E <- function(h2, E_ratio) {
    var_E_m <- get_varE_m(h2)
    var_E_f <- get_varE_f(E_ratio, var_E_m)
    E <- NULL
    for (i in 1:length(sex)) {
        if (sex[i] == 0) {
            E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_f[i]))
        } else {
            E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_m[i]))
        }
    }
    return(E)
}

gwas_mash <- function(E) {
    ### GWAS ###
    pheno <- colSums(Beta*snps) + E
    f_model <- lm(pheno[f_id] ~ snp[f_id,] + E[f_id])
    m_model <- lm(pheno[m_id] ~ snp[m_id,] + E[m_id])
    f_gwas <- summary(f_model)$coefficients[1:snp_num+1,1:2]
    m_gwas <- summary(m_model)$coefficients[1:snp_num+1,1:2]
    mash_BETA <-  matrix(c(f_gwas[,1],m_gwas[,1]), nrow=poly_num, ncol=2)
    mash_SE <-  matrix(c(f_gwas[,2],m_gwas[,2]), nrow=poly_num, ncol=2)
    colnames(mash_BETA) <- c('female','male') ; colnames(mash_SE) <- c('female','male')

    ### mash ###
    data <- mash_set_data(mash_BETA, mash_SE)
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

h2_list <- c(0.5, 0.05)
E_ratio_list <- c(1,1.5,5)
for (E_ratio in E_ratio_list) {
    for (h2 in h2_list) {
        E <- get_E(h2,E_ratio)
        gwas_mash(E)
    }
}



