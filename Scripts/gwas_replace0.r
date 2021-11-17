#!/usr/bin/env Rscript
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno

wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno)
setwd(wd)

remove_zero <- function(sex) {
    file_name <- paste0(sex,"_all.",pheno,".glm.linear")
    df <- read.table(file_name, sep="\t", head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses = c(rep("integer",2), rep("character", 6), "integer", rep("numeric", 4  )))

    zeros <- df[df$P ==0,]
    if (nrow(zeros) > 0) {
        write.table(zeros, file=paste0(pheno,"_glm_zeros_",sex,".txt"),sep="\t",row.names=FALSE)
    }
    min_p <- min(df[df$P != 0,]$P)
    df$P[df$P==0] <- min_p
    write.table(df, file=paste0(sex,"_all_no0.",pheno,".glm.linear"),quote=FALSE)
}

remove_zero("male")
remove_zero("female")
remove_zero("both_sex")