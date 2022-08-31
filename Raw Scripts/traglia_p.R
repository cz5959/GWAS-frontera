#!/usr/bin/env Rscript

# TEST ENVIRONMENTAL VARIANCE IN MASH
library(optparse, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")

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

# sex p
get_p <- function(beta, se) {
  t <- beta / se
  p <- 2*pt(abs(t), 300000-2, lower.tail = FALSE)
  return(p)
}
mash_BETA$p_f <- get_p(mash_BETA$female, mash_SE$female)
mash_BETA$p_m <- get_p(mash_BETA$male, mash_SE$male) 

# traglia sex het
z <- (mash_BETA$female - mash_BETA$male) / sqrt( mash_SE$male^2 + mash_SE$female^2)
p <- 2*pnorm(abs(z), 0, 1, lower.tail=FALSE)
p_id <- which(p<5e-8)

# sex heterogenous snps
BETA_sexhet <- mash_BETA[p_id,]
print(paste("number sex het snps: ", nrow(BETA_sexhet)))

# sex heterogenous snps - nonsig
nonsig <- which((BETA_sexhet$p_f < 5e-8 & BETA_sexhet$p_m < 5e-8))
print(paste("nonsig: ", length(one_sex)))

# effect in only one sex; one sex sig one sex not
one_sex <- which((BETA_sexhet$p_f < 5e-8 & BETA_sexhet$p_m > 5e-8) | (BETA_sexhet$p_f > 5e-8 & BETA_sexhet$p_m < 5e-8))
print(paste("effect in only one sex: ", length(one_sex)))
# effect in opposite directions; both sex sig or not sig, diff sign
opp <- which((BETA_sexhet$female < 0 & BETA_sexhet$male > 0) | (BETA_sexhet$female > 0 & BETA_sexhet$male < 0))
opp <- opp[!(opp %in% one_sex)]
print(paste("effect in opposite directions: ", length(opp)))
print(paste("     opposite, nonsig: ", length(opp[(opp %in% nonsig)])))

# differences in magnitude of effect; both sex sig or not sig, same sign
same <- which((BETA_sexhet$female > 0 & BETA_sexhet$male > 0) | (BETA_sexhet$female < 0 & BETA_sexhet$male < 0))
same <- same[!(same %in% one_sex)]
print(paste("effect in same direction: ", length(same)))
print(paste("     same, nonsig: ", length(same[(same %in% nonsig)])))