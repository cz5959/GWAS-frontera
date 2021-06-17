#!/usr/bin/env Rscript

setwd("/scratch1/08005/cz5959/Association_Height_50")

library("qqman")

results_as <- read.table("linear_results_all_chrom.height.glm.linear",sep="\t",head=FALSE,col.names=c("#CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))

jpeg("manhattan_height.jpeg")

manhattan(results_as,chr="#CHROM",bp="POS",p="P",snp="SNP", main = "Manhattan plot: Height")
dev.off()  




