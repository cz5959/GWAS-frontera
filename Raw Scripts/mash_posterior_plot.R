#!/usr/bin/env Rscript

libLoc <- "/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/"
library("optparse", lib.loc=libLoc)
library("crayon",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)
library("tidyr",lib.loc=libLoc)
library("withr",lib.loc=libLoc)
library("ggplot2",lib.loc=libLoc)
library("grid",lib.loc=libLoc)
library("gridExtra",lib.loc=libLoc)

option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)

wd <- paste0("/scratch/08005/cz5959/GWAS_Results/",pheno,"/mash")
setwd(wd)

# sample indexes
set.seed(1)
# load mash file for gwas estimates and LD blocks
load(file= paste0(pheno,"_mash.RData"))
# sample LD blocks
random <- numeric(0)
for (i in unique(LD_groups$group)) {
    sample_subset <- LD_groups[LD_groups$group == i, 'index']
    random[[(length(random) + 1)]] <- sample(sample_subset,1)
}      
print(paste0("random subset length: ", length(random)))
return(random)

# GWAS F-M estimate regression
gwas_df <- data$Bhat; rm(data); rm(LD_groups)
gwas_df <- as.data.frame(gwas_df)
gwas_model <- lm(female ~ male, gwas_df)
gwas_B <- summary(gwas_model)$coefficients[2]; gwas_y_i <- summary(gwas_model)$coefficients[1]
gwas_r2 <- summary(gwas_model)$r.squared
gwas_df <- gwas_df[random,]        

# PM
# mash posterior mean F-M regression
pm_df <- read.csv(paste0(pheno,"_mash_pm.txt"), sep="\t", colClasses = rep("numeric",2))
pm_model <- lm(female ~ male, pm_df)
pm_B <- summary(pm_model)$coefficients[2]; pm_y_i <- summary(pm_model)$coefficients[1]
pm_r2 <- summary(pm_model)$r.squared
pm_df <- pm_df[random,]  

# mash mixture proportions line; slope foe beta_f v. beta_m for SNPs sampling across matrices
weight <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t"); weight <- weight[c(2:101)]
weight <- rowMeans(weight)
setwd("/scratch/08005/cz5959/mash_posterior_plots")
ratio <- read.csv("matrice_names.txt", sep="\t", colClasses=c(rep("NULL",4), "numeric")); ratio <- ratio$f_to_m_magnitude
null_weight <- sum(weight[c(1,3,4)]); null_weight <- 1/(1-null_weight)
weight <- weight[-c(1,3,4)]*null_weight ; ratio <- ratio[-c(1,3,4)]  # remove null
alpha <- sum(weight * ratio)

# print slopes
print( paste0(pheno, "\n raw slope: ", gwas_B, "\n pm slope: ", pm_B))

# plot
setwd("/scratch/08005/cz5959/mash_posterior_plots")
xmax <- max( max(gwas_df$male), max(pm_df$male) )
xmin <- min( min(gwas_df$male), min(pm_df$male) )
ymax <- max( max(gwas_df$female), max(pm_df$female) )
ymin <- min( min(gwas_df$female), min(pm_df$female) )

pdf(file=paste0(pheno,"_pm_gwas.pdf"), width=8, height=4)

pm_label <- paste0("slope = ", round(pm_B,3))
pm_label2 <- bquote(R^2 == .(format(pm_r2, digits = 3)))
mix_label <- paste0("alpha (mash) = ", round(alpha,3))
p1 <- ggplot(pm_df, aes(x=male, y=female)) +        # posterior mean plot
    geom_point(size=0.5, alpha = 0.5) +
    geom_abline(slope = pm_B, intercept = pm_y_i, size=0.8, alpha = 0.5, color= "black") +
    geom_abline(slope = alpha, intercept = 0, size=0.8, alpha = 0.5, color="#2b62d9") +
    theme_classic() +
    labs(x="Male Posterior Mean", y="Female Posterior Mean", title="Posterior effect estimates") +
    scale_x_continuous(limits=c(xmin,xmax)) + scale_y_continous(limits=c(ymin,ymax)) +
    theme(plot.title = element_text(size=13), axis.title = element_text(size=11), axis.text = element_text(size=9))
pm_grob <- grobTree(textGrob(pm_label, x=0.1, y=0.9, hjust=0, gp=gpar(col="black", fontsize=9)))
pm_grob2 <- grobTree(textGrob(pm_label2, x=0.1, y=0.84, hjust=0, gp=gpar(col="black", fontsize=9)))
mix_grob <- grobTree(textGrob(mix_label, x=0.1, y=0.78, hjust=0, gp=gpar(col="#2b62d9", fontsize=9)))
p1 <- p1 + annotation_custom(pm_grob) + annotation_custom(pm_grob2) + annotation_custom(mix_grob)

gwas_label <- paste0("slope = ", round(gwas_B,3))
gwas_label2 <- bquote(R^2 == .(format(gwas_r2, digits = 3)))
p2 <- ggplot(gwas_df, aes(x=male, y=female)) +      # raw gwas estimate plot
    geom_point(size=0.5, alpha = 0.5) +
    geom_abline(slope = gwas_B, intercept = gwas_y_i, size=0.8, alpha = 0.5, color= "black") +
    theme_classic() +
    labs(x="Male Raw Effect Estimate", y="Female Raw Effect Estimate", title="Raw GWAS effect estimates") +
    scale_x_continuous(limits=c(xmin,xmax)) + scale_y_continous(limits=c(ymin,ymax)) +
    theme(plot.title = element_text(size=13), axis.title = element_text(size=11), axis.text = element_text(size=9)) 
gwas_grob <- grobTree(textGrob(gwas_label, x=0.1, y=0.9, hjust=0, gp=gpar(col="black", fontsize=9)))
gwas_grob2 <- grobTree(textGrob(gwas_label2, x=0.1, y=0.84, hjust=0, gp=gpar(col="black", fontsize=9)))
p2 <- p2 + annotation_custom(gwas_grob) + annotation_custom(gwas_grob2)

p <- grid.arrange(p2, p1, nrow = 1)
print(p)

dev.off()