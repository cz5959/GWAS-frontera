setwd("~/Research/GWAS-frontera/mash")

pheno <- "height"
df <- read.csv(paste0(pheno,"_mash_partitioned.txt"), sep="\t")