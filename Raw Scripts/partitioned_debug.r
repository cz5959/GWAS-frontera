#!/usr/bin/env Rscript
libLoc <- "/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/"
library(optparse, lib.loc=libLoc)
library("crayon",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)

## examine overlap of cell type groups ##
setwd("/scratch1/08005/cz5959/LD_practice/Partitioned/celltype_to_snp")
snp_df <- read.csv("celltype_all.txt", sep="\t", colClasses=c(rep("integer",12)))
ttl_snps <- nrow(snp_df)        # 9607691

snp_df2 <- snp_df %>%            
    select(-c(X.CHROM,POS)) %>% 
    filter(!if_all(starts_with("Group"), ~ . == 0)) %>% # remove rows with all zeros
    mutate_if(is.numeric, ~1 * (. != 0)) # convert all non-zero to 1
group_snps <- nrow(snp_df)        # 4231304

snp_sum <- rowSums(snp_df)
celltype_sum <- colSums(snp_df)

write.table(snp_sum, file="celltype_sums.txt",row.names=FALSE)
write.table(celltype_sum, file="celltype_colsums.txt",sep="\t",row.names=FALSE)

## col SUM posterior weights to see if match mixture weights ##
pheno <- 'height'
setwd(paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash"))
df <- read.csv(paste0(pheno,"_mash_weights.txt"),sep="\t",nrow=1)
num <- length(df)
df <- read.csv(paste0(pheno,"_mash_weights.txt"),sep="\t",nrows=9607691,colClasses = c(rep("NULL",0), rep("numeric",5), rep("NULL",(num-5))))
df_sums <- colSums(df)
mixdf <- read.csv(paste0(pheno,"mixprop_100_all.txt",sep="\t"))






