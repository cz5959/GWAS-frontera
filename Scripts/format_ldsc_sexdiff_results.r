#!/usr/bin/env Rscript

library(optparse, lib.loc="/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/")
# parameters
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-s","--sex"), type="character", default="sexdiff",help="sex",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)
sex <- opt$sex; print(sex)

# load results and h2 file
setwd(paste0("/scratch/08005/cz5959/LD_practice/",pheno,"/celltypes"))
df <- read.csv(paste0(pheno,"_",sex,"_celltype_enrichment.txt"), sep="\t", head = FALSE,
col.name=c("Type","Enrichment","Enrichment.SE"))
h2 <- read.csv(paste0(pheno,"_",sex,"_celltype_h2.txt"), sep=" ", head = FALSE,
col.name=c("a","b","c","d","h2","h2.se"),
colClasses=c(rep("NULL",4),"numeric","character"))
h2$h2.se <- as.numeric(gsub("[()]", "", h2$h2.se))
# cell types
celltypes <- c("Skeletal Muscle", "Adrenal Pancreas", "Cardiovascular", "CNS", "Connective Bone",
"GI", "Hematopoietic", "Kidney", "Liver", "Other")
df$Type = celltypes

# merge h2 with enrichment df
df <- cbind(df,h2)

write.table(df, paste0(pheno,"_",sex,"_celltype_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)