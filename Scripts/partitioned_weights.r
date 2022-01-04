#!/usr/bin/env Rscript
libLoc <- "/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/"
library(optparse, lib.loc=libLoc)
library("crayon",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)

option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)

# load cell type list
setwd("/scratch1/08005/cz5959/LD_practice/Partitioned/celltype_to_snp")
snp_df <- read.csv("celltype_all.txt", sep="\t", colClasses=c(rep("integer",12)))

# load chromosome and positions 
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno)
setwd(wd)
pos_df <- read.csv(paste0("female_all.",pheno,".glm.linear"), sep="\t", colClasses = c(rep("integer",2), rep("NULL", 11)))
colnames(pos_df)=c("CHROM","POS")

# get list of covariance matrices for phenotype
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd)
df <- read.csv(paste0(pheno,"_mash_weights.txt"),sep="\t",nrow=1)
num <- length(df)

get_weight <- function(df) {
    # combine position with df
    df <- data.frame( cbind(pos_df, df) )   
    # combine weights with cell groups
    df <- df %>%
        merge(snp_df, by.x=c("CHROM","POS"), by.y=c("X.CHROM","POS")) %>%
        select(-c(CHROM,POS))
    # group by cell type, get mean of each weight
    final_grouped <- data.frame(matrix(ncol=6,nrow=0))
    for (i in 1:10) {
        grouped_df <- df %>%
            group_by(across(5+i)) %>%
            summarise(across(c(1:5), list(sum)))
        grouped_df <- data.frame(grouped_df)
        final_grouped[nrow(final_grouped)+1,] <- grouped_df[2,] 
    }
    colnames(final_grouped) <- c("group",colnames(df[1:5]))
    return(final_grouped[,2:6])
}

# iterate over sections of columns in mash weights files
final_df <- data.frame(matrix(ncol=0,nrow=10))
for (i in 1:ceiling(num/5)) {
    print(i)
    beg = (i-1)*5
    end = num - (5 + beg)
    if (i == ceiling(num/5) & num %% 5 > 0) {
        df <- read.csv(paste0(pheno,"_mash_weights.txt"),sep="\t",nrows=9607691,colClasses = c(rep("NULL",beg), rep("numeric",num%%5), rep("NULL",0)))
    } else {
        df <- read.csv(paste0(pheno,"_mash_weights.txt"),sep="\t",nrows=9607691,colClasses = c(rep("NULL",beg), rep("numeric",5), rep("NULL",end)))
    }
    final_df <- cbind(final_df, get_weight(df))
    save(final_df, file=paste0(pheno,i,"_mash_partitioned.RData"))
}

write.table(final_df, file=paste0(pheno,"_mash_partitioned.txt"), sep="\t", row.names=FALSE)


####
perl -pe 's/(\d*\.\d*)/sprintf("%.9f",$1+0.5)/ge' testosterone_mash_weights.txt > testosterone_mash_weights2.txt
