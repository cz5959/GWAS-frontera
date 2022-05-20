#!/usr/bin/env Rscript
libLoc <- "/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/"
library(optparse, lib.loc=libLoc)
library("crayon",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)

option_list = list(
    make_option(c("-n","--num"), type="character", default=NULL,help="cell type",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
num <- opt$num; print(num)

# load snps (height)
#setwd("/scratch1/08005/cz5959/GWAS_Results/height/")
#snp_df <- read.csv("female_all.height.glm.linear", sep="\t", colClasses=c(rep("integer",2),rep("NULL",11)))
#snp_df <- snp_df %>%
#    mutate(group = 0) %>%
#    arrange(X.CHROM)

# open cell type file
get_celltype <- function(num) {
    setwd("/scratch1/08005/cz5959/LD_practice/Partitioned/1000G_Phase3_cell_type_groups")
    cell_df <- read.csv(paste0(num,".bed"), sep="\t", colClasses=c("character",rep("integer",2)), col.names=c("chr","start","end"))
    cell_df <- cell_df %>%
        mutate(chr=substr(chr,4,nchar(chr)+1)) %>%
        filter(!(chr %in% c("X","Y","M")))
    setwd("/scratch1/08005/cz5959/LD_practice/Partitioned/celltype_to_snp")
    return(cell_df)
}

#
setwd("/scratch1/08005/cz5959/LD_practice/Partitioned/celltype_to_snp")
#### load wip snps (type=1)
load("celltype_to_snp.RData")
snp_df <- select(snp_df, -group)

####
print(paste("celltype: ",num))
snp_df$group <- 0
cell_df <- get_celltype(num)
for (i in 1:22) {           # split by chromosome
    cell <- cell_df[cell_df$chr == i,]
    id = 1  # index of cell type groups
    for (j in which(snp_df$X.CHROM == i)) {       # snps of specified chromosome
        if (id > nrow(cell)) { break }      # check if reached the end of cell group
        pos <- snp_df$POS[j]
        while (pos > cell$start[id] & id <= nrow(cell)) {       # check if greater than start
            if (pos <= cell$end[id]) {      # assign cell group if within range
                snp_df$group[j] <- num
                break
            } else {        # move to next position range
                id <- id + 1
            }
        }
    }
}
#save(snp_df, file=paste0(num,"_celltype_to_snp_wip.RData"))

save(snp_df, file=paste0(num,"_celltype_to_snp.RData"))
write.table(snp_df, file=paste0(num,"_celltype_to_snp.txt"), sep="\t", row.names=FALSE)
    

