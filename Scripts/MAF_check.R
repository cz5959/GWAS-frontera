args <- commandArgs(trailingOnly = TRUE)
n <- args[1]
title <- sprintf("MAF distribution Chr%s",n)
maf_freq <- read.table("MAF_check.afreq", header =TRUE, as.is=T)
pdf(sprintf("MAF_distribution_%s.pdf",n))
hist(maf_freq[,5],main = title, xlab = "ALT_FREQS")
dev.off()


