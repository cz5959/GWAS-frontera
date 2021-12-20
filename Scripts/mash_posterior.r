#!/usr/bin/env Rscript

# load libraries
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-m","--mode"), type="character", default="additive",help="mash for PGS?",metavar="character"),
    make_option(c("-s","--set"), type="character", default="1",help="set number",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)
set <- opt$set; print(set)

# add PGS suffix if using PGS method 
if (opt$mode == "additive") {
    suffix <- "" ; wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")

} else { 
    suffix <- "_pgs" ; wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash_",set)
} 
setwd(wd)
load(file= paste0(pheno,"_mash",suffix,".RData"))
load(file= paste0(pheno,"_mash_100g",suffix,".RData"))

# adjust table for mixture weights
weight_col <- function(df) {
    colnames(df) <- gsub("^(.*)[.].*", "\\1",colnames(df))
    df <- t(rowsum(t(df), group = colnames(df)))
    return(df)
}

# posterior summaries for all
header <- c("female", "male")
pm_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(pm_all) <- header
psd_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(psd_all) <- header
lfsr_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(lfsr_all) <- header

# split calculation into chunks of 30k
interval <- 30000
num <- floor(nrow(data$Bhat) / interval)
for (i in 0:num) {
    print(paste0("progress: ", i,"/",num))
    start <- (i*interval) + 1
    end <- (i+1)*interval
    if (end > nrow(data$Bhat)) {
        end <- nrow(data$Bhat)
    } 
    datasub=mash_set_data(data$Bhat[start:end,], data$Shat[start:end,])
    msub = mash(datasub, g=g_ave, fixg=TRUE)    
    if (i == 0) {
        cov_names <- colnames(msub$posterior_weights)
        weights_all <- data.frame(matrix(ncol = length(cov_names), nrow = 0)) ; colnames(weights_all) <- cov_names
    }
    pm = get_pm(msub)
    psd = get_psd(msub)
    lfsr = get_lfsr(msub)
    pm_all <- rbind(pm_all, pm)
    psd_all <- rbind(psd_all, psd)
    lfsr_all <- rbind(lfsr_all, lfsr)
    #save(pm_all,psd_all,lfsr_all,file=paste0(pheno,"_mash_posterior",suffix,".RData"))
    weights <- weight_col(msub$posterior_weights)
    weights_all <- rbind(weights_all, msub$posterior_weights)
    #save(pm_all,psd_all,lfsr_all,weights_all, file=paste0(pheno,"_mash_posterior",suffix,".RData"))
}

# write posterior estimates to table
write.table(pm_all, file=paste0(pheno,"_mash_pm",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(psd_all, file=paste0(pheno,"_mash_psd",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(lfsr_all, file=paste0(pheno,"_mash_lfsr",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(weights_all, file=paste0(pheno,"_mash_weights",suffix,".txt"), sep="\t", row.names=FALSE)