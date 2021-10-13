#!/usr/bin/env Rscript

### mash lfsr to pseudo-pvalues ###

# argument parser
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load mash lfsr and pm results
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd)
lfsr_df <- read.table(paste0(pheno,"_mash_lfsr_pgs.txt"),sep="\t",head=TRUE,colClasses=c(rep("numeric",2)))
pm_df <- read.table(paste0(pheno,"_mash_pm_pgs.txt"),sep="\t",head=TRUE,colClasses=c(rep("numeric",2)))

wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/PGS")
setwd(wd)
get_pseudo_p <- function(sex) {
    # load gwas p-values, 2 columns: ID and P
    file_name <- paste0(sex,"_train.",pheno,".glm.logistic")
    m_gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses=c(rep("NULL",2), "character", rep("NULL", 2), "character", rep("NULL",6), "numeric")) 

    # sort p-values
    sorted_p <- gwas_df[order(gwas_df$P),3]
    # concatenate ids, A1, lfsr and pm
    if (sex == "female") {
        mash_ids <- data.frame(gwas_df$ID, gwas_df$A1, pm_df$female, lfsr_df$female)
    } else {
        mash_ids <- data.frame(gwas_df$ID, gwas_df$A1, pm_df$male, lfsr_df$male)
    }
    
    colnames(mash_ids) <- c("ID", "A1", "pm", "lfsr")
    # sort by lfsr 
    mash_ids <- mash_ids[order(mash_ids$lfsr),]
    # concatenate sorted ids with pseudo pvalues
    mash_ids$P <- sorted_p

    # write to text file
    write.table(mash_ids, file=paste0(sex,"_pseudoP_pgs.",pheno,".txt"), sep="\t", row.names=FALSE, quote=FALSE)

    return(mash_ids)
}

# plot lfsr against p-values
for (sex in c("female","male")) {
    df <- get_pseudo_p(sex)
    pdf( paste0("lfsr-p_",pheno,"_",sex,".pdf") )
    df_small <- df[c(rep(FALSE,99),TRUE),]  # get every 100th SNP
    plot(lfsr ~ P, data=df_small, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), xlab="p-values", ylab="lfsr", pch=".")

    pdf( paste0("lfsr-p_",pheno,"_",sex,"_5e-2.pdf") )
    plot(lfsr ~ P, data=df_small, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), xlab="p-values", ylab="lfsr", pch=".",xlim=c(0,0.05),ylim=c(0,0.05))
    png( paste0("lfsr-p_",pheno,"_",sex,"_1e-3.png") )
    plot(lfsr ~ P, data=df, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), xlab="p-values", ylab="lfsr", pch=".",xlim=c(0,1e-3),ylim=c(0,1e-3))
    png( paste0("lfsr-p_",pheno,"_",sex,"_1e-5.png") )
    plot(lfsr ~ P, data=df, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), xlab="p-values", ylab="lfsr", pch=".",xlim=c(0,1e-5),ylim=c(0,1e-5))
}

df <- read.table(paste0(sex,"_pseudoP_pgs_small.",pheno,".txt"),sep="\t",header=TRUE,colClasses=c(rep("NULL",3),rep("numeric",2)))
