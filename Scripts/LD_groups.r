#!/usr/bin/env Rscript

# 2 minutes

# load in variants 
setwd("/scratch1/08005/cz5959/QC")
df <- read.table("ukb_imp_all_v3_11.pvar",sep="\t",head=FALSE, col.names=c("CHROM","POS","ID","REF","ALT"), colClasses = c(rep("integer",2), rep("NULL", 3)))
df$index <- 1:nrow(df)

# Pickrell LD blocks
setwd("/scratch1/08005/cz5959/LD_practice/LD_scores")
load(file="LD_blocks_EUR.RData")

# LD_scores
LD_blocks$index <- 1:nrow(LD_blocks)
pos_groups <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pos_groups) <- c("CHROM", "POS","index","group")
for (i in 1:22) {
    chr <- df[df$CHROM == i,]
    chr$group <- 0
    ld <- LD_blocks[LD_blocks$chr == i,]
    for (j in 1:nrow(ld)) {
        chr$group <- ifelse( chr$POS >= ld$start[j] & chr$POS < ld$stop[j], ld$index[j], chr$group)
    }
    pos_groups <- rbind(pos_groups, chr)
}
save(pos_groups,file="LD_groups.RData")