# TEST ENVIRONMENTAL VARIANCE IN MASH

library(ashr)
library(mashr)
library(ggplot2)
library(reshape2)
setwd("~/Research/GWAS-frontera/mash/simulation")

### 350k ## if statement to remove monomorphism ######################heat plots

# create sample of individuals -- random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

### SNP ###
#### 170,000 pickrell LD blocks into mash
# do the 100 causal first, then do gwas sites (can ignore the causal)

# set seed
set.seed(1)

# get genotype at site
get_genotypes <- function(snp_freq) {
  allele1 <- rbinom(n=n, size=1, prob=snp_freq)
  allele2 <- rbinom(n=n, size=1, prob=snp_freq)
  genotypes <- allele1 + allele2
  return(genotypes)
}

# i causal SNPs
snp_num <- 100 ; snp_num_name <- '1e2'
snp_freqs <- read.csv(paste0("maf_sample_",snp_num_name,".txt"), colClasses = 'numeric') ; snp_freqs <- snp_freqs$x

genotype_matrix_i <- matrix(nrow=n, ncol=length(snp_freqs))
monomorphs <- NULL
for (i in 1:length(snp_freqs)) {
  snp_count <- get_genotypes(snp_freqs[i])
  if (length(unique(snp_count)) ==1) {       # check if snp is monomorphic
    monomorphs[i] <- i
    next
  }
  #genotype_matrix_i[,i] <- snp_count
  genotype_matrix_i <- cbind(genotype_matrix_i, snp_count)
  if (i %% 1000 == 0){
    print(i)
  }
}

######################################################################################

# Beta (effect size)
Beta <- rnorm(snp_num, mean=0, sd=1)

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G <- (Beta^2) * snp_freqs * (1-snp_freqs)   # genetic variance
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

# environment and heritability combination
for (h2 in c(0.5, 0.05)) {
  for (E_ratio in c(1, 1.5, 5)) {
    E <- get_environment(h2, E_ratio)
    pheno <- (Beta * rowSums(genotype_matrix_i)) + E
    assign(paste0("E_",h2,"_",E_ratio), E)
    assign(paste0("pheno_",h2,"_",E_ratio), pheno)
  }
}

# GWAS ## TODo: more snps per gwas
GWAS <- function(genotype, pheno, E) {
  f_model <- lm(pheno[f_id] ~ genotype[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ genotype[m_id] + E[m_id])
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

################################
# k snps, LD blocks
gwas_snps <- read.csv("maf_sample_1700e2.txt", colClasses = 'numeric') ; gwas_snps <- gwas_snps$x
h2 = 0.5; E_ratio = 1
E <- get_environment(h2, E_ratio)
pheno <- (Beta * rowSums(genotype_matrix_i)) + E

genotype_matrix_k <- NULL
for (i in 1:length(gwas_snps)) {
  snp <- get_genotypes(gwas_snps[i])    # get allele count for all individuals
  genotype_matrix_k <- cbind(genotype_matrix_k, snp)
  if (i %% 1000 == 0){
    print(i)
  }
}

############# load genotype matrix file 10000 at a time

mash_BETA <- NULL
mash_SE <- NULL
for (i in 1:170000) {
  gwas <- GWAS(genotype_matrix_k[,1], E)
  mash_BETA <- rbind(mash_BETA, data.frame(female= gwas[1,1], male= gwas[2,1]))
  mash_SE <- rbind(mash_SE, data.frame(female= gwas[1,2], male= gwas[2,2]))
  if (i %% 1000 == 0){
    print(i)
  }
}

mash_func <- function() {
  for (h2 in c(0.5, 0.05)) {
    for (E_ratio in c(1, 1.5, 5)) {
      BETA <- mash_BETA[mash_BETA$h2 == h2 & mash_BETA$E_ratio == E_ratio, c(3,4)]
      SE <- mash_SE[mash_SE$h2 == h2 & mash_SE$E_ratio == E_ratio, c(3,4)]
      BETA <- as.matrix(BETA) ; SE <- as.matrix(SE)
      data <- mash_set_data(BETA, SE, zero_Shat_reset = .Machine$double.eps)
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

# check BETA and SE
if (length(monomorphs) == length(snpids)) {
  print("All SNPs are monomorphic :(")
} else {
  # run the mash function
  mash_func()
} 


#############################################
mix_names <- names(U.c)
mix_names <- c("identity",mix_names)

mix <- mix_names
for (h2 in c(0.5,0.05)) {
  for (envr in c(1, 1.5, 5)) {
    new_mix <- read.csv(paste0("100_",h2,"_",envr,".txt"), sep="\t")
    colnames(new_mix) <- paste0("100_",h2,"_",envr)
    mix <- cbind(mix, new_mix)
  }
}
mix_long <- melt(mix, id.vars="mix")

ggplot(data=mix_long, aes(x=mix, y=value)) +
  geom_point(size=2) +
  facet_wrap(~variable, ncol=6) +
  coord_flip()


