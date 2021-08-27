#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# load arguments
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
print(pheno)
wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",pheno,sep="")
setwd(wd)
load(file= paste(pheno,"_mash.RData",sep=""))

# random subset: sample once from every LD block
set.seed(1)
random <- numeric(0)
for (i in unique(LD_groups$group)) {
    sample_subset <- LD_groups[LD_groups$group == i, 'index']
    random[[(length(random) + 1)]] <- sample(sample_subset,1)
}
print(paste("random subset ", length(random), sep=""))

#correlation structure
data.temp = mash_set_data(data$Bhat[random,],data$Shat[random,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)

# set up canoncial covar matrices and add hypothesis
U.c = cov_canonical(data.random)
# add hypothesis
corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
effect = c(1.5,2,3)
for (c in corr) {
    for (e in effect) {
        U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
        U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
    }
}
U.c[['opp_het_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
U.c[['opp_het_2']] <- matrix(c(1,-0.5,-0.5,1),2,2)
U.c[['opp_het_3']] <- matrix(c(1,-0.75,-0.75,1),2,2)
U.c[['equal_opp']] <- matrix(c(1,-1,-1,1),2,2)
names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")

# fit mash model 
m = mash(data.random, Ulist = U.c, outputlevel = 1)

# posterior summaries for all
pm_all <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(pm_all) <- c("female", "male")
psd_all <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(psd_all) <- c("female", "male")
lsfr_all <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(lsfr_all) <- c("female", "male")
interval <- 30000
num <- floor(nrow(data$Bhat) / interval)
for (i in 0:num) {
    start <- (i*interval) + 1
    end <- (i+1)*interval
    if (end > nrow(data$Bhat)) {
        end <- nrow(data$Bhat)
    } 
    datasub=mash_set_data(data$Bhat[start:end,], data$Shat[start:end,])
    msub = mash(datasub, g=get_fitted_g(m), fixg=TRUE)
    pm = get_pm(msub)
    psd = get_psd(msub)
    lfsr = get_lfsr(msub)
    pm_all <- rbind(pm_all, pm)
    psd_all <- rbind(psd_all, psd)
    lsfr_all <- rbind(lsfr_all, lfsr)
}

write.table(pm_all, file=paste(pheno,"_mash_pm.txt",sep=""), sep="\t", row.names=FALSE)
write.table(psd_all, file=paste(pheno,"_mash_psd.txt",sep=""), sep="\t", row.names=FALSE)
write.table(lsfr_all, file=paste(pheno,"_mash_lfsr.txt",sep=""), sep="\t", row.names=FALSE)

save(m_strong, m, file= paste(pheno,"_mash_results.RData",sep=""))