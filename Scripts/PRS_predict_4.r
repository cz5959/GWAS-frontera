#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
type <- args[2]
print(pheno)

p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)

# load phenotype and covariate file
setwd("/scratch1/08005/cz5959/Phenotypes")
pheno_file <- paste0("pheno_",pheno,".txt")
phenotype <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("integer","integer","numeric"))
covariates <- read.table("covariates.txt", header=FALSE, colClasses=c(rep("integer",2), rep("numeric", 12), rep("NULL", 3)))
colnames(covariates)=c("FID","IID",paste0("PC",1:10),"sex","birthyear")
# merge phenotype and covariate file
pheno_covar <- merge(phenotype, covariates, by=c("FID", "IID"))

PRS_prediction <- function(sex) {  
    # test set
    test_ids <- read.table(paste0(pheno,"_",sex,"_testIIDs.txt"), sep="\t")
    pheno_covar <- pheno_covar[pheno_covar$IID %in% test_ids,]
    # null model
    null.model <- lm(pheno ~ ., data=pheno_covar[,!colnames(pheno_covar) %in% c("FID","IID")])
    null.r2 <- summary(null.model)$r.squared
    null.coeff <- summary(null.model)$coefficients

    prs.result <- NULL
    prs.residuals <- NULL
    for(i in p.threshold){
        # go through each p-value threshold
        prs <- read.table(paste0(sex,"_",type,"_",pheno,".",i,".profile"), header=T)
        # slice into test set ids       ########### should already be only test ids
        prs <- prs[prs$IID %in% test_ids,]
        # Merge prs with phenotype matrix - only take FID, IID and PRS 
        pheno.prs <- merge(pheno_covar, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))

        # linear regression
        model <- lm(pheno~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
        model.r2 <- summary(model)$r.squared        # model R2
        prs.r2 <- model.r2-null.r2                  # incremental R2
        # other summary results                                           
        prs.coef <- summary(model)$coefficients["SCORE",]
        prs.residuals <- summary(model)$residuals
        prs.beta <- as.numeric(prs.coef[1])
        prs.se <- as.numeric(prs.coef[2])
        prs.p <- as.numeric(prs.coef[4])
        # store results
        prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=model.r2, incR2=prs.r2, P=prs.p, BETA=prs.beta, SE=prs.se))
        prs.residuals <- rbind(prs.residuals, summary(model)$residuals)
    }
    # Best threshold based on R2 is:
    print(paste0("Best result for ", sex, ", ",type))
    print(prs.result[which.max(prs.result$R2),])
    # return results
    prs_summary <- list("results" = prs.result, "residuals" = prs.residuals)
    return(prs_summary)
}

female_prs.summary <- PRS_prediction("female")
male_prs.summary <- PRS_prediction("male")

# save results as Rdata
save(female_prs.summary, male_prs.summary, file= paste0(pheno,"_",type,"_prs_results.RData"))