#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# load arguments
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
mode <- args[2]
print(pheno)
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd)

if (is.na(mode)) { suffix <- "" } else { suffix <- "_prs" } 
load(file= paste0(pheno,"_mash",suffix".RData"))
load(file= paste0(pheno,"_mash_100g",suffix".RData"))

# posterior summaries for all
header <- c("female", "male")
pm_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(pm_all) <- header
psd_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(psd_all) <- header
lsfr_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(lsfr_all) <- header

interval <- 30000
num <- floor(nrow(data$Bhat) / interval)
for (i in 0:num) {
    start <- (i*interval) + 1
    end <- (i+1)*interval
    if (end > nrow(data$Bhat)) {
        end <- nrow(data$Bhat)
    } 
    datasub=mash_set_data(data$Bhat[start:end,], data$Shat[start:end,])
    msub = mash(datasub, g=g_ave, fixg=TRUE)
    pm = get_pm(msub)
    psd = get_psd(msub)
    lfsr = get_lfsr(msub)
    pm_all <- rbind(pm_all, pm)
    psd_all <- rbind(psd_all, psd)
    lsfr_all <- rbind(lsfr_all, lfsr)
}

write.table(pm_all, file=paste0(pheno,"_mash_pm",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(psd_all, file=paste0(pheno,"_mash_psd",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(lsfr_all, file=paste0(pheno,"_mash_lfsr",suffix,".txt"), sep="\t", row.names=FALSE)