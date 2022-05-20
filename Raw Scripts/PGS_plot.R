
pheno <- "height"


load(file=paste0(pheno,"_additive_pgs_results.RData"))
female_add <- female_pgs.summary
male_add <- male_pgs.summary
load(file=paste0(pheno,"_",type,"_both_sex_pgs_results.RData"))
female_add_both <- female_pgs.summary
male_add_both <- male_pgs.summary
load(file=paste0(pheno,"_mash_pgs_results.RData"))
female_mash <- female_pgs.summary
male_mash <- male_pgs.summary

pdf( paste0(pheno,"_",type,"_female_pgs_residuals.pdf") )
actual <- female_add$fitted[i,] - female_add$residuals[i,]
plot(actual, female_add$fitted[i,], main=paste("PGS: Additive and mash: female, ",pheno,sep=" "), 
    col='red', xlab="y", ylab="y_hat", pch=".")
points(actual, female_mash$fitted[i,], col='green')
points(add, female_add_both$fitted[i,], col='blue')
legend(x="topright",legend=c("Add. Same","Add. Both","mash Same"),col=c("red","blue","green"))

pdf( paste0(pheno,"_",type,"_male_pgs_residuals.pdf") )
actual <- male_add$fitted[i,] - male_add$residuals[i,]
plot(actual, male_add$fitted[i,], main=paste("PGS: Additive and mash: male, ",pheno,sep=" "), 
    col='red',xlab="y", ylab="y_hat", pch=".")
points(actual, male_mash$fitted[i,],col='green')
points(add, male_add_both$fitted[i,],col='blue')
legend(x="topright",legend=c("Add. Same","Add. Both","mash Same"),col=c("red","blue","green"))