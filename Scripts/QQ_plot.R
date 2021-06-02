#!/usr/bin/env Rscript

# set working directory
#setwd("/scratch1/08005/cz5959/Association_Height_50")
setwd("/scratch1/08005/cz5959/Neale_Lab/Standing_Height")

# install package qqman
install.packages("qqman", repos="https://cran.microsoft.com/", lib="~")
library("qqman")

# read results and convert P to numeric
#results_log <- read.table("linear_results_all_chrom.height.glm.linear",sep="\t",head=FALSE,col.names=c("#CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","TSTAT","P"))
#P <- as.numeric(results_log$P)

#Neale
results_log <- read.table("50_raw.gwas.imputed_v3.both_sexes.tsv", sep="\t", head = FALSE, col.names=c("variant","minor_allele","minor_AF","low_confidence_variant","n_complete_sample","AC","ytx","beta","se","tstat","pval"))
P <- as.numeric(results_log$pval)
jpeg("QQ-Plot_Linear_Height_Neale.jpeg")
qq(P, main = "Q-Q plot of Neale GWAS p-values for Height")

#create qq plot
#jpeg("QQ-Plot_Linear_Height.jpeg")
#qq(P, main = "Q-Q plot of GWAS p-values for Height")

dev.off()


