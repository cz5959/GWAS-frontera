#!/usr/bin/env Rscript

# arguments and set working directory
setwd("/scratch1/08005/cz5959/Phenotypes")
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
print(pheno)

################ get samples from merged QC file instead
# IID dataframes
pheno_file <- paste0("pheno_",pheno,".txt")
df_pheno <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("NULL","integer","numeric"))
df_sex <- read.table("sex_ids.txt", sep="\t", head=FALSE, col.names=c("IID","sex"))        
df_sex_check <- read.table("diffsex_ids.txt", sep="\t", head=FALSE, col.names=c("IID"))     #############
df_sex_aneuploidy <- read.table("sexaneuploidy_ids.txt", sep="\t", head=FALSE, col.names=c("IID"))      #################
df_relatives <- ####################
df_wb = read.table("whitebritishIDs.txt", sep=" ",head=FALSE,col.names=c("FID","IID"), colClasses=c("NULL","integer"))

# merge sex, and pheno
df <- merge(df_pheno, df_sex, by="IID")

# keep white british
df <- df[df$IID %in% df_wb$IID,]

# remove relatives ######################

# remove diff sex, and sex aneuploidy ids
df <- df[df$IID %in% df_sex_check$IID,]
df <- df[df$IID %in% df_sex_aneuploidy$IID,]

# split into female=0 and male=1
female_df <- df[df$sex == 0,]
male_df <- df[df$sex == 1]

# sample 25k
set.seed(1)
f_test <- sample(female_df,25000)
m_test <- sample(male_df,25000)

wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/PRS")
setwd(wd)
write.table(f_test$IID, file=paste0(pheno,"_female_testIIDs.txt"), sep="\t", row.names=FALSE)
#write.table(f_test, file=paste0(pheno,"_female_test.txt"), sep="\t", row.names=FALSE)
write.table(m_test$IID, file=paste0(pheno,"_male_testIIDs.txt"), sep="\t", row.names=FALSE)
#write.table(m_test, file=paste0(pheno,"_male_test.txt"), sep="\t", row.names=FALSE)









