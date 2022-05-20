#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# arguments and set working directory
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
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

strong_f <- female_df[female_df$ID %in% female_clump$SNP,]
strong_f <- strong_f[order(strong_f$P),]
strong_f <- strong_f[!duplicated(strong_f$ID),]
strong_m <- male_df[male_df$ID %in% male_clump$SNP,]
strong_m <- strong_m[order(strong_m$P),]
strong_m <- strong_m[!duplicated(strong_m$ID), ]


# new ID column CHROM:POS:REF:ALT:ID
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
# select strong signals
strong_f <- which(female_df$P < 5e-8)
strong_m <- which(male_df$P < 5e-8)
strong <- unique(c(strong_f,strong_m))
print(paste("sig snps > 5e-8: ", length(strong), sep=""))
strong_f <- rownames(female_df)
strong_m <- rownames(male_df)
strong <- unique(c(strong_f,strong_m))
print(paste("clumped snps ", length(strong), sep=""))
#m.1by1 = mash_1by1(data)
#strong = get_significant_results(m.1by1,0.05)
# get random subset
random = sample( 1:nrow(BETA), (nrow(BETA)/10) )
#correlation structure
data.temp = mash_set_data(BETA[random,],SE[random,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)
data.strong = mash_set_data(data$Bhat[strong,],data$Shat[strong,], V=Vhat)
# obtain initial data-driven covar matrices (PCA) and apply extreme deconvolution
U.pca = cov_pca(data,2,subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)
U.pca = cov_pca(data.strong,2)
U.ed = cov_ed(data.strong, U.pca)
# set up canoncial covar matrices
U.c = cov_canonical(data)
U.c = cov_canonical(data.random)
# fit mash model 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
# compute posterior summmaries
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE) # for strong tests
# run mash using both canonical and datadriven covariances
m = mash(data, c(U.c,U.ed))

#### ANALYSIS ####
# extract posterior summaries
posteriors <- cbind(get_lfsr(m2), get_pm(m2), get_psd(m2))
colnames(posteriors) <- c("F_LFSR", "M_LFSR", "F_PM", "M_PM", "F_PSD", "M_PSD")
write.table(posteriors, file = paste(pheno,"_posteriors_strong.txt",sep=""), sep = "\t")
# sharing
sharing <- get_pairwise_sharing(m2) # sharing by magnitude
write.table(sharing, file = paste(pheno,"sharing_strong.txt",sep=""), sep = "\t")
# measure of fit, log-likelihood
print(paste("log-likelihood:",get_loglik(m2),sep="")
# mixture proportions bar plot
mixture_prop <- get_estimated_pi(m2)
write.table(mixture_prop, file=paste(pheno,"mixprop_strong.txt",sep=""), sep-"\t")
png( paste(pheno,"_mixture_bar_strong.png",sep=""))
par(mar=c(7,4,4,2)+.1)
barplot(mixture_prop,las=2,main= paste("Mixture Proportions - ",pheno,sep=""), ylab="Proportion")
pdf( paste(pheno,"_mixture_bar_strong.pdf",sep=""))
barplot(mixture_prop,las=2,,main= paste("Mixture Proportions - ",pheno,sep=""), ylab="Proportion")) #las=2 perpindicular labels


dev.off()

