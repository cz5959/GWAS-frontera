#!/usr/bin/env Rscript

# load libraries
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-m","--mode"), type="character", default="additive",help="mash for PGS?",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno

# load LD scores and set working directory
setwd("/scratch1/08005/cz5959/LD_practice/LD_scores")
load(file="LD_groups.RData")
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno)
setwd(wd)

setup_sex_df <- function(sex) {
    # load results - summstats
    if (opt$mode == "additive") {
        file_name <- paste0(sex,"_all.",pheno,".glm.linear")
    } else {
        file_name <- paste0(wd,"/PGS/",sex,"_train.",pheno,".glm.linear")
    }
    
    gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric")) 

    # add new ID column CHROM:POS:REF:ALT:ID
    gwas_df$VAR <- paste(gwas_df$CHROM, gwas_df$POS, gwas_df$REF, gwas_df$ALT, gwas_df$ID, sep=":")
    gwas_df$index <- seq.int(nrow(gwas_df))

    return(gwas_df)
}

## binary
setup_sex_df <- function(sex) {
    # load results - summstats
    file_name <- paste0(wd,"/PGS/",sex,"_train.",pheno,".glm.logistic")
    
    gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric")) 

    # add new ID column CHROM:POS:REF:ALT:ID
    gwas_df$VAR <- paste(gwas_df$CHROM, gwas_df$POS, gwas_df$REF, gwas_df$ALT, gwas_df$ID, sep=":")
    gwas_df$index <- seq.int(nrow(gwas_df))

    return(gwas_df)
}

####

female_df <- setup_sex_df("female")
male_df <- setup_sex_df("male")

####

# LD groups merge with p-values; pos_groups from loaded LD scores
LD_groups <- merge(pos_groups, female_df[c("index","P")], by="index")
LD_groups <- merge(LD_groups, male_df[c("index","P")], by="index")

####

# create matrix
conditions <- c("female", "male")
r <- nrow(female_df)
BETA <- matrix(c(female_df$BETA, male_df$BETA), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))
SE <- matrix(c(female_df$SE, male_df$SE), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))

# create mash data object
data = mash_set_data(BETA, SE)

# save variables into Rdata
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd) 
if (opt$mode != "additive") {
    save(data, LD_groups, file= paste0(pheno,"_mash_pgs.RData"))
} else {
    summstat_pos <- female_df$POS
    save(summstat_pos, LD_groups, data, file= paste0(pheno,"_mash.RData"))
}

################################# mismatch in variants in logistic regression
female_df <- female_df[female_df$ID %in% male_df$ID,]
male_df <- male_df[male_df$ID %in% female_df$ID,]
female_df$ld_match <- paste(female_df$CHROM,female_df$POS,sep=":")
male_df$ld_match <- paste(male_df$CHROM,male_df$POS,sep=":")
pos_groups$ld_match <- paste(pos_groups$CHROM,pos_groups$POS,sep=":")
female_df$index <- seq.int(nrow(female_df))
male_df$index <- seq.int(nrow(male_df))
LD_groups <- merge(pos_groups, female_df[c("ld_match","index","P")], by="ld_match")
LD_groups <- merge(LD_groups, male_df[c("ld_match","index","P")], by="ld_match")
LD_groups <- LD_groups[c(2,3,5,7,8,9)]