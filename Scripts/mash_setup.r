#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# arguments and set working directory
setwd("/scratch1/08005/cz5959/LD_practice/LD_scores")
load(file="LD_groups.RData")
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
mode <- args[2]
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno)
setwd(wd)

setup_df <- function(sex) {
    # load results - summstats
    if (is.na(mode)) {
        file_name <- paste0(sex,"_all.",pheno,".glm.linear")
    }
    else {
        file_name <- paste0(wd,"/PRS/",sex,"_train.",pheno,".glm.linear")
    }
    
    gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric")) 

    #load results - clumped
    clumped_file_name <- paste0(sex,"_",pheno,"_tab.clumped")
    clumped_df <- read.table(clumped_file_name,sep="\t",head=TRUE)

    # add new ID column CHROM:POS:REF:ALT:ID
    gwas_df$VAR <- paste(gwas_df$CHROM, gwas_df$POS, gwas_df$REF, gwas_df$ALT, gwas_df$ID, sep=":")
    gwas_df$index <- seq.int(nrow(gwas_df))

    # LD groups merge with p-values
    LD_groups <- merge(pos_groups, gwas_df[c("index","P")], by="index")

    return(gwas_df)
}

female_df <- setup_df("female")
male <- setup_df("female")

# LD groups merge with p-values
LD_groups <- merge(pos_groups, female_df[c("index","P")], by="index")
LD_groups <- merge(LD_groups, male_df[c("index","P")], by="index")

# create matrix
conditions <- c("female", "male")
r <- nrow(female_df)
BETA <- matrix(c(female_df$BETA, male_df$BETA), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))
SE <- matrix(c(female_df$SE, male_df$SE), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))

# create mash data object
data = mash_set_data(BETA, SE)

wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd) 
if (! is.na(mode)) {
    save(data, LD_groups, file= paste0(pheno,"_mash_prs.RData"))
}
else {
    # strong subset index list
    strong_f <- female_df[female_df$ID %in% female_clump$SNP,]
    strong_f <- strong_f[strong_f$P < 5e-8,'index']
    strong_m <- male_df[male_df$ID %in% male_clump$SNP,]
    strong_m <- strong_m[strong_m$P < 5e-8,'index']
    strong <- unique(c(strong_f,strong_m))

    # save Rdata
    summstat_pos <- df$POS

    save(summstat_pos, LD_groups, data, strong, file= paste0(pheno,"_mash.RData"))
}