#!/usr/bin/env Rscript

wd <- paste("/scratch1/08005/cz5959/LD_practice/LD_scores")
setwd(wd)
load(file="LD_blocks_EUR.RData")

file <- "both_sex_all.height.glm.linear"
df <- read.table(file,sep="\t",head=FALSE, 
col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric"))

# LD_scores
summstat_pos <- df[, c('CHROM','POS','index')]
LD_blocks$index <- 1:nrow(LD_blocks)
pos_groups <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(pos_groups) <- c("CHROM", "POS","index","group")
for (i in 1:22) {
    chr <- summstat_pos[summstat_pos$CHROM == i,]
    chr$group <- 0
    ld <- LD_blocks[LD_blocks$chr == i,]
    for (j in 1:nrow(ld)) {
        chr$group <- ifelse( chr$POS >= ld$start[j] & chr$POS < ld$stop[j], ld$index[j], chr$group)
    }
    pos_groups <- rbind(pos_groups, chr)
}
save(pos_groups,file="LD_groups.RData")