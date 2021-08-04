#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# arguments and set working directory
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]

wd <- paste("/scratch1/08005/cz5959/LD_practice/LD_scores")
setwd(wd)
load(file="LD_blocks_EUR.RData")

wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",pheno,sep="")
setwd(wd)

# load results - summstats
female_file <- paste("female_all.",pheno,".glm.linear",sep="")
male_file <- paste("male_all.",pheno,".glm.linear",sep="")
female_df <- read.table(female_file,sep="\t",head=FALSE, 
col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric"))
male_df <- read.table(male_file,sep="\t",head=FALSE,
col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"),
colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric"))

#load results - clumped
female_file <- paste("female_",pheno,"_tab.clumped",sep="")
male_file <- paste("male_",pheno,"_tab.clumped",sep="")
female_clump <- read.table(female_file,sep="\t",head=TRUE)
male_clump <- read.table(male_file,sep="\t",head=TRUE)

# add new ID column CHROM:POS:REF:ALT:ID
female_df$VAR <- paste(female_df$CHROM, female_df$POS, female_df$REF, female_df$ALT, female_df$ID, sep=":")
male_df$VAR <- paste(male_df$CHROM, male_df$POS, male_df$REF, male_df$ALT, male_df$ID, sep=":")

# create matrix
conditions <- c("female", "male")
r <- nrow(female_df)
BETA <- matrix(c(female_df$BETA, male_df$BETA), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))
SE <- matrix(c(female_df$SE, male_df$SE), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))

#### MASH ####
# read in data
data = mash_set_data(BETA, SE)

# strong subset index list
strong_f <- female_df[female_df$ID %in% female_clump$SNP,]
strong_f <- which(strong_f$P < 5e-8)
strong_m <- male_df[male_df$ID %in% male_clump$SNP,]
strong_m <- which(strong_m$P < 5e-8)
strong <- unique(c(strong_f,strong_m))

# LD_scores
summstat_pos <- female_df$POS
# summstat_pos <- female_df[, c('CHROM','POS')]
# LD_blocks$index <- 1:nrow(LD_blocks)

# pos_groups <- data.frame(matrix(ncol = 3, nrow = 0))
# colnames(pos_groups) <- c("CHROM", "POS", "group")
# for (i in 1:22) {
#     chr <- summstat_pos[summstat_pos$CHROM == i,]
#     chr$group <- 0
#     ld <- LD_blocks[LD_blocks$chr == i,]
#     for (j in 1:nrow(ld)) {
#         chr$group <- ifelse( chr$POS >= ld$start[j] & chr$POS < ld$stop[j], ld$index[j], chr$group)
#     }
#     pos_groups <- rbind(pos_groups, chr)
# }
# save(pos_groups,file=,"LD_groups.RData")

# save Rdata
save(summstat_pos, data, strong, file= paste(pheno,"_mash.RData",sep=""))