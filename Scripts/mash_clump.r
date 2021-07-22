#!/usr/bin/env Rscript

# arguments and set working directory
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",pheno,sep="")
setwd(wd)
print(pheno)
load(file= paste(pheno,"_mash.RData",sep=""))
#### MASH ####

# strong subset
print(paste("strong subset ", length(strong), sep=""))
# random subset from sampling
random = sample( 1:nrow(BETA), (nrow(BETA)/50) )
# random subset from partitioned sampling
female_df <- female_df[order(female_df$CHROM,female_df$POS),]
start = female_df
print(paste("random subset ", length(random), sep=""))

#correlation structure
data.temp = mash_set_data(BETA[random,],SE[random,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)
data.strong = mash_set_data(data$Bhat[strong,],data$Shat[strong,], V=Vhat)

# obtain initial data-driven covar matrices (PCA) and apply extreme deconvolution
U.pca = cov_pca(data.strong,2)
print(str(U.pca))
U.ed = cov_ed(data.strong, U.pca)

# set up canoncial covar matrices and add hypothesis
U.c = cov_canonical(data.random)
# add hypothesis
U.c[['f_1']] <- f_1 <- matrix(c(1.25,0,0,1),2,2)
U.c[['f_2']] <- f_2 <- matrix(c(2,0,0,1),2,2)
U.c[['f_3']] <- f_3 <- matrix(c(4,0,0,1),2,2)
U.c[['f_corr_1']] <- f_corr_1 <- matrix(c(2.25,1.5,1.5,1),2,2)
U.c[['f_corr_2']] <- f_corr_2 <- matrix(c(4,2,2,1),2,2)
U.c[['f_corr_3']] <- f_corr_3 <- matrix(c(9,3,3,1),2,2)
U.c[['m_1']] <- m_1 <- matrix(c(1,0,0,1.25),2,2)
U.c[['m_2']] <- m_2 <- matrix(c(1,0,0,2),2,2)
U.c[['m_3']] <- m_3 <- matrix(c(1,0,0,4),2,2)
U.c[['m_corr_1']] <- m_corr_1 <- matrix(c(1,1.5,1.5,2.25),2,2)
U.c[['m_corr_2']] <- m_corr_2 <- matrix(c(1,2,2,4),2,2)
U.c[['m_corr_3']] <- m_corr_3 <- matrix(c(1,3,3,9),2,2)
U.c[['f_het_1']] <- f_het_1 <- matrix(c(2.25,0.75,0.75,1),2,2)
U.c[['f_het_2']] <- f_het_2 <- matrix(c(4,1,1,1),2,2)
U.c[['f_het_3']] <- f_het_2 <- matrix(c(9,1.5,1.5,1),2,2)
U.c[['m_het_1']] <- m_het_1 <- matrix(c(1,0.75,0.75,2.25),2,2)
U.c[['m_het_2']] <- m_het_1 <- matrix(c(1,1,1,4),2,2)
U.c[['m_het_3']] <- m_het_1 <- matrix(c(1,1.5,1.5,9),2,2)
U.c[['f_opp']] <- f_opp <- matrix(c(2.25,-1.5,-1.5,1),2,2)
U.c[['m_opp']] <- m_opp <- matrix(c(1,-1.5,-1.5,2.25),2,2)
U.c[['het_opp_1']] <- het_opp_1 <- matrix(c(1,-0.25,-0.25,1),2,2)
U.c[['het_opp_2']] <- het_opp_2 <- matrix(c(1,-0.5,-0.5,1),2,2)
U.c[['equal_opp']] <- equal_opp <- matrix(c(1,-1,-1,1),2,2)

# fit mash model 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

#### ANALYSIS ####
# mixture proportions bar plot
mixture_prop <- get_estimated_pi(m)
write.table(mixture_prop, file=paste(pheno,"mixprop.txt",sep=""), sep="\t")
png( paste(pheno,"_mixture_bar.png",sep=""), width = 1200, height = 480)
par(mar=c(7,4,4,2)+.1)
barplot(mixture_prop,las=2,main= paste("Mixture Proportions - ",pheno,sep=""), ylab="Proportion")

dev.off()