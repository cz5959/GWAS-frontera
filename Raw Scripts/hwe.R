args <- commandArgs(trailingOnly = TRUE)
n <- args[1]

title1 <- sprintf("Histogram HWE Chr%s",n)
title2 <- sprintf("Histogram HWE: strongly deviating SNPs only; Chr%s",n)

hwe<-read.table (file="plink2.hardy", header=TRUE)
pdf(sprintf("histhwe_%s.pdf",n))
hist(hwe[,10],main=title1)
dev.off()

hwe_zoom<-read.table (file="plink2zoomhwe.hardy", header=TRUE)
pdf(sprintf("histhwe_below_threshold_%s.pdf",n))
hist(hwe_zoom[,10],main=title2)
dev.off()
