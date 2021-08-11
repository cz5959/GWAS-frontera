#!/usr/bin/env Rscript

# load libraries
library(ashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library(mashr, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# arguments and set working directory
args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
print(pheno)
wd <- paste("/scratch1/08005/cz5959/LD_practice/LD_scores")
setwd(wd)
load(file="LD_groups.RData")
wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",pheno,sep="")
setwd(wd)
load(file= paste(pheno,"_mash.RData",sep=""))

#### MASH ####

# strong subset
print(paste("strong subset ", length(strong), sep=""))

random_subset <- function(seed=1) {
    set.seed(seed)

    # METHOD 1: random sampling
    # random = sample( 1:length(summstat_pos), (length(summstat_pos)/50) )
    # from partitioned sampling, sample once per 250000kb

    # METHOD 2: sample every ~250,000 kb
    num <- floor((max(summstat_pos) - min(summstat_pos)) / 250000)
    interval <- floor(length(summstat_pos) / num)
    random <- numeric(0)
    for (i in 0:num) {
        start <- i*interval
        if (start+interval > length(summstat_pos)) {
            random[[(length(random) + 1)]] <- sample(start:length(summstat_pos),1)
        } else {
            random[[(length(random) + 1)]] <- sample(start:start+interval,1)
        }
    }

    # METHOD 3: sample once from every LD block
    random <- numeric(0)
    for (i in unique(pos_groups$group)) {
        sample_subset <- which(pos_groups$group == i)
        random[[(length(random) + 1)]] <- sample(sample_subset,1)
    }
    return(random)

    print(paste("random subset ", length(random), sep=""))
    return(random)
}

fit_mash <- function(random) {
    #correlation structure
    data.temp = mash_set_data(data$Bhat[random,],data$Shat[random,])
    Vhat = estimate_null_correlation_simple(data.temp)
    rm(data.temp)
    data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)
    data.strong = mash_set_data(data$Bhat[strong,],data$Shat[strong,], V=Vhat)

    # obtain initial data-driven covar matrices (PCA) and apply extreme deconvolution
    U.pca = cov_pca(data.strong,2)
    U.ed = cov_ed(data.strong, U.pca)
    print(str(U.pca))
    print(str(U.ed))

    # set up canoncial covar matrices and add hypothesis
    U.c = cov_canonical(data.random)
    # add hypothesis
    corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
    effect = c(1.5,2,3)
    for (c in corr) {
        for (e in effect) {
            U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
            U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
        }
    }

    U.c[['opp_het_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
    U.c[['opp_het_2']] <- matrix(c(1,-0.5,-0.5,1),2,2)
    U.c[['opp_het_3']] <- matrix(c(1,-0.75,-0.75,1),2,2)
    U.c[['equal_opp']] <- matrix(c(1,-1,-1,1),2,2)

    # fit mash model 
    m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
    # mixture proportions bar plot
    mixture_prop <- get_estimated_pi(m)
    return(mixture_prop)
}

for (i in 1:20) {
    random <- random_subset(i)
    ## sample w/o replacement
    pos_groups <- pos_groups[-c(random),]
    
    mix <- fit_mash(random)
    if (i==1) {
        mixture <- matrix(names(mix), ncol=1)
    }
    mixture <- cbind(mixture, mix)
}
colnames(mixture) <- cbind("names",paste("mix_",1:20,sep=""))
write.table(mixture, file=paste(pheno,"mixprop.txt",sep=""), sep="\t", row.names=FALSE)
#png( paste(pheno,"_mixture_bar.png",sep=""), width = 1200, height = 480)
#par(mar=c(7,4,4,2)+.1)
#barplot(mixture_prop,las=2,main= paste("Mixture Proportions - ",pheno,sep=""), ylab="Proportion")

#dev.off()