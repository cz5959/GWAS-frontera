#!/bin/sh

TARGET=$SCRATCH/1000G/EUR

# QC target data
plink --bfile $TARGET --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out ${TARGET}.QC

# prune to remove highly correlated SNPs
plink --bfile $TARGET --keep ${TARGET}.QC.fam --extract ${TARGET}.QC.snplist --indep-pairwise 200 50 0.25 --out ${TARGET}.QC

# compute heterozygosity rates
plink --bfile $TARGET --extract ${TARGET}.QC.prune.in --keep ${TARGET}.QC.fam --het --out ${TARGET}.QC

### R ### manual
dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R

# sex check
plink --bfile $TARGET --extract ${TARGET}.QC.prune.in --keep ${TARGET}.valid.sample --check-sex --out ${TARGET}.QC

### R ### manual
# Read in file
valid <- read.table("EUR.valid.sample", header=T)
dat <- read.table("EUR.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
q() # exit R

# relatedness 
plink --bfile $TARGET --extract ${TARGET}.QC.prune.in --keep ${TARGET}.QC.valid --rel-cutoff 0.125 --out ${TARGET}.QC

# final QC data set
plink --bfile $TARGET --make-bed --keep ${TARGET}.QC.rel.id --out ${TARGET}.QC --extract ${TARGET}.QC.snplist