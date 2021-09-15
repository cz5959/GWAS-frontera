#!/usr/bin/env Rscript

### mash lfsr to pseudo-pvalues ###
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
print(pheno)

# load mash lfsr and pm results
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd)
lfsr_df <- read.table(paste0(pheno,"_mash_lsfr_prs.txt"),sep="\t")
pm_df <- read.table(paste0(pheno,"_mash_pm_prs.txt"),sep="\t")

wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/PRS")
setwd(wd)
get_pseudo_p <- function(sex) {
    # load gwas p-values, 2 columns: ID and P
    file_name <- paste0(sex,"_train.",pheno,".glm.linear")
    gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses=c(rep("NULL",2), "character", rep("NULL", 2), "character", rep("NULL",6), "numeric")) 

    # sort p-values
    sorted_p <- gwas_df[order(gwas_df$P),3]
    # concatenate ids, A1, lfsr and pm
    mash_ids <- cbind(gwas_df$ID, gwas_df$A1, pm_df$sex, lfsr_df$sex)
    colnames(mash_ids) <- c("ID", "A1", "pm", "lfsr")
    # sort by lfsr 
    mash_ids <- mash_ids[order(mash_ids$lfsr),]
    # concatenate sorted ids with pseudo pvalues
    pseudo_p <- cbind(mash_ids, sorted_p)

    # write to text file
    write.table(pseudo_p, file=paste0(sex,"_pseudoP.",pheno,".txt"), sep="\t", row.names=FALSE)

    return(pseudo_p)
}

# plot lfsr against p-values
for (sex in c("female","male")) {
    df <- get_pseudo_p(sex)
    pdf( paste0("lfsr-p_",pheno,"_",sex,".pdf") )
    plot(lfsr ~ P, data=df, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), xlab="p-values", ylab="lfsr", pch=".")
}
