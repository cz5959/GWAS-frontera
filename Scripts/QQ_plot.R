#!/usr/bin/env Rscript

# set working directory
setwd("/scratch1/08005/cz5959/Association_Height_50")

# install package qqman
install.packages("qqman", repos="https://cran.microsoft.com/", lib="~")
library("qqman")

# create qq plot
results_log <- read.table("linear_results_all_chrom.height.glm.linear",head=TRUE)
jpeg("QQ-Plot_Linear_Height.jpeg")
qq(results_log$P, main = "Q-Q plot of GWAS p-values for Height")
dev.off()


