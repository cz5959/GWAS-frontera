#!/usr/bin/env Rscript

# arguments and names
args <- commandArgs(trailingOnly=TRUE)
wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",args[1],sep="")
both_file <- paste("both_sex_all.",args[1],".glm.linear",sep="")
female_file <- paste("female_all.",args[1],".glm.linear",sep="")
male_file <- paste("male_all.",args[1],".glm.linear",sep="")
# set working directory
setwd(wd)

# qqman package
#install.packages("qqman", repos="https://cran.microsoft.com/", lib="~")
library("qqman", lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# create QQ plot: read results, convert P to numeric, create jpeg
results_log <- read.table(both_file,sep="\t",head=FALSE,col.names=c("#CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
P <- as.numeric(results_log$P)
jpeg( paste("qqplot_",args[1],"_both.jpeg",sep="") )
qq(P, main = paste("QQ Plot of GWAS p-values for ",toupper(args[1]),": Both Sex", sep="") )

results_log <- read.table(female_file,sep="\t",head=FALSE,col.names=c("#CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
P <- as.numeric(results_log$P)
jpeg( paste("qqplot_",args[1],"_female.jpeg",sep="") )
qq(P, main = paste("QQ Plot of GWAS p-values for ",toupper(args[1]),": Female", sep="") )

results_log <- read.table(male_file,sep="\t",head=FALSE,col.names=c("#CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
P <- as.numeric(results_log$P)
jpeg( paste("qqplot_",args[1],"_male.jpeg",sep="") )
qq(P, main = paste("QQ Plot of GWAS p-values for ",toupper(args[1]),": Male", sep="") )

#Neale
#setwd("/scratch1/08005/cz5959/Neale_Lab/Standing_Height")
#results_log <- read.table("50_raw.gwas.imputed_v3.both_sexes.tsv", sep="\t", head = FALSE, col.names=c("variant","minor_allele","minor_AF","low_confidence_variant","n_complete_sample","AC","ytx","beta","se","tstat","pval"))
#P <- as.numeric(results_log$pval)
#jpeg("QQ-Plot_Linear_Height_Neale.jpeg")
#qq(P, main = "Q-Q plot of Neale GWAS p-values for Height")

dev.off()


