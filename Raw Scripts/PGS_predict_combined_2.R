#!/usr/bin/env Rscript

# argument parser
library(optparse, lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-t","--type"), type="character", default=NULL,help="additive or mash",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)
type <- opt$type; print(type)

p.threshold <- c('1', '0.01', '1e-5', '1e-8')

# load phenotype and covariate file
setwd("/scratch1/08005/cz5959/Phenotypes")
pheno_file <- paste0("pheno_",pheno,".txt")
phenotype <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("integer","integer","numeric"))
covariates <- read.table("covariates.txt", header=FALSE, colClasses=c(rep("integer",2), rep("numeric", 12), rep("NULL", 3)))
colnames(covariates)=c("FID","IID",paste0("PC",1:10),"sex","birthyear")
# merge phenotype and covariate file
pheno_covar <- merge(phenotype, covariates, by=c("FID", "IID"))

# test set
wd <- paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno,"/PGS")
setwd(wd)
f_ids <- read.table(paste0(pheno,"_female_testIIDs.txt"), sep="\t", col.names=c("FID","IID"))
m_ids <- read.table(paste0(pheno,"_male_testIIDs.txt"), sep="\t", col.names=c("FID","IID"))
test_ids <- rbind(f_ids,m_ids)
pheno_covar <- pheno_covar[pheno_covar$IID %in% test_ids$IID,]

# null model
null.model <- lm(paste0(pheno," ~ ."), data=pheno_covar[,!colnames(pheno_covar) %in% c("FID","IID")])
print(summary(null.model))
null.r2 <- summary(null.model)$r.squared
null.coeff <- summary(null.model)$coefficients

pgs.result <- NULL
pgs.residuals <- NULL
pgs.fitted <- NULL
for(i in p.threshold){
    print(i)
    # go through each p-value threshold
    pgs <- read.table(paste0("both_sex_",type,"_",pheno,".",i,".profile"), header=T)
    pgs <- pgs[pgs$IID %in% test_ids$IID,]
    pgs.snps <- mean(pgs$CNT)/2
    # Merge pgs with phenotype matrix - only take FID, IID and pgs 
    pheno.pgs <- merge(pheno_covar, pgs[,c("FID","IID", "SCORE")], by=c("FID","IID"))

    # linear regression
    model <- lm(paste0(pheno," ~ ."), data=pheno.pgs[,!colnames(pheno.pgs)%in%c("FID","IID")])
    model.r2 <- summary(model)$r.squared        # model R2
    pgs.r2 <- model.r2-null.r2                  # incremental R2
    # other summary results                                           
    pgs.coef <- summary(model)$coefficients["SCORE",]
    pgs.beta <- as.numeric(pgs.coef[1])
    pgs.se <- as.numeric(pgs.coef[2])
    pgs.p <- as.numeric(pgs.coef[4])
    # store results
    pgs.result <- rbind(pgs.result, data.frame(Threshold=i, R2=model.r2, incR2=pgs.r2, P=pgs.p, BETA=pgs.beta, SE=pgs.se))
    pgs.residuals <- rbind(pgs.residuals, summary(model)$residuals)
    pgs.fitted <- rbind(pgs.fitted, model$fitted.values)
}

# Best threshold based on R2 is:
print(paste0("Best result for combined - both sex, ", type))
print(pgs.result[which.max(pgs.result$R2),])
# return results
combined_pgs.summary <- list("results" = pgs.result, "residuals" = pgs.residuals, "fitted" = pgs.fitted)