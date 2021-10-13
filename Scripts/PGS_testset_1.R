#!/usr/bin/env Rscript

# argument parser
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)

set.seed(1)

# set working directory
setwd("/scratch1/08005/cz5959/Phenotypes")

################ get samples from merged QC file instead
# IID dataframes
pheno_file <- paste0("pheno_",pheno,".txt")
df_ids <- read.table("QC_ids.txt", sep="\t", head=FALSE, col.names=c("IID")) 
df_pheno <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("NULL","integer","numeric"))
df_sex <- read.table("sex_ids.txt", sep="\t", head=FALSE, col.names=c("IID","sex"))        

# merge sex, and pheno
df <- merge(df_pheno, df_sex, by="IID")

# keep QC ids
df <- df[df$IID %in% df_ids$IID,]

sets=5
for (i in 1:sets) {
    # split into female=0 and male=1
    female_df <- df[df$sex == 0,]
    male_df <- df[df$sex == 1,]

    # sample 25k and convert to dataframe
    f_test <- sample(female_df$IID,25000)
    f_test <- data.frame(f_test,f_test)
    m_test <- sample(male_df$IID,25000)
    m_test <- data.frame(m_test,m_test)

    # save as file
    wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/PGS_",i)
    setwd(wd)
    write.table(f_test, file=paste0(pheno,"_female_testIIDs.txt"), sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(m_test, file=paste0(pheno,"_male_testIIDs.txt"), sep="\t", row.names=FALSE, col.names=FALSE)

    # remove those test ids from df
    df <- df[! df$IID %in% f_test$f_test,]
    df <- df[! df$IID %in% m_test$m_test,]
}



