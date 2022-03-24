#!/usr/bin/env Rscript

library(optparse, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
# parameters
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)

# load male and female gwas file
setwd(paste0("/scratch/08005/cz5959/GWAS_Results/",pheno))

gwas_df <- read.csv(paste0("male_all.",pheno,".glm.linear"), sep="\t", head=TRUE,
    colClasses = c(rep("NULL",2), "character", rep("NULL", 2), rep("character", 2), "NULL", "integer", rep("numeric", 2), rep("NULL", 2)) ) 
f_df <- read.csv(paste0("female_all.",pheno,".glm.linear"), sep="\t", head=TRUE,
    colClasses = c(rep("NULL",2), "character", rep("NULL", 2), rep("character", 2), "NULL", "integer", rep("numeric", 2), rep("NULL", 2)) ) 

# check if gwas is of different lengths
if (nrow(gwas_df) != nrow(f_df)) {
    print('ERROR: different gwas lengths!')
    if (nrow(gwas_df) > nrow(f_df)) {
        print( gwas_df[!gwas_df$ID %in% f_df$ID,] )
        gwas_df <- gwas_df[gwas_df$ID %in% f_df$ID,]
    } else {
        print( f_df[!f_df$ID %in% gwas_df$ID,] )
        f_df <- f_df[f_df$ID %in% gwas_df$ID,]
    }
}

# beta = male_beta - female_beta
# se = sqrt ( male_se^2 + female_se^2)
# z = beta / se
# p = pnorm(zscore, mean=0, sd=1)
gwas_df$BETA <- gwas_df$BETA - f_df$BETA 
gwas_df$SE <- sqrt( (gwas_df$SE^2) + (f_df$SE^2) )
gwas_df$Z <- gwas_df$BETA / gwas_df$SE
gwas_df$P <- pnorm(gwas_df$Z, 0, 1)

write.table(gwas_df, paste0("sex_diff.",pheno,".glm.linear"), sep="\t", row.names=FALSE, quote=FALSE)






