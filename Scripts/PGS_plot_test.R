###########################################


height_f=c(0.2056,0.1590,0.2142)
height_m=c(0.2075,0.1441,0.2028)
bmi_f=c(0.0604,0.0346,0.0582)
bmi_m=c(0.0546,0.0346,0.0566)
testosterone_f=c(0.0005,0.0082,0.0076)
testosterone_m=c(0.0652,0.0636,0.0627)
waist_hip_bmi_f=c(0.0425,0.0516,0.0490)
waist_hip_bmi_m=c(0.0161,0.0135,0.0166)
df <- cbind(height_f,height_m,bmi_f,bmi_m,testosterone_f,testosterone_m,waist_hip_bmi_f,waist_hip_bmi_m)
types<-c('Add-both',"Add-same","mash-same")
rownames(df) <- types
head(df)
colors=c('azure1','azure2','azure3')
barplot(df,beside=T,col=colors,ylim=c(0,0.26),legend=TRUE,args.legend=(list(x="top",ncol=3)))

#########################
sex<-"female";pheno<-"height"
df <- read.table(paste0(sex,"_pseudoP_pgs_small.",pheno,".txt"),
                 sep="\t",header=TRUE,colClasses=c(rep("NULL",3),rep("numeric",2)))
head(df)
limit<-1
plot(lfsr ~ P, data=df, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), 
     xlab="p-values", ylab="lfsr", pch=".",xlim=c(0,limit),ylim=c(0,limit))

plot(lfsr ~ P, data=df, main=paste("lfsr vs p-values -",sex,pheno,sep=" "), 
     xlab="p-values", ylab="lfsr", pch=".",log="x",ylim=c(0,1e-8))

plot(lfsr ~ (log10(P)), data=df, main=paste("lfsr vs log10(p-values) -",sex,pheno,sep=" "), 
     xlab="log10(p-values)", ylab="lfsr", pch=".",ylim=c(0,1e-5))


#########################
pheno<-"height"
load(file=paste0(pheno,"_additive_pgs_results.RData"))
female_add <- female_pgs.summary
male_add <- male_pgs.summary
load(file=paste0(pheno,"_additive_both_sex_pgs_results.RData"))
female_add_both <- female_pgs.summary
male_add_both <- male_pgs.summary
load(file=paste0(pheno,"_mash_pgs_results.RData"))
female_mash <- female_pgs.summary
male_mash <- male_pgs.summary

# female additive - same
i=2
male_add_both$results
actual <- female_add$fitted[i,] - female_add$residuals[i,]
plot(actual, female_add$fitted[i,], main=paste("PGS: Additive and mash: female, ",pheno,sep=" "), 
     col='red', xlab="y", ylab="y_hat", pch=".")
# female additive - both
i=1
points(actual, female_add_both$fitted[i,], col='blue',pch=".")
# female mash
i=1
points(actual, female_mash$fitted[i,], col='green',pch=".")

#height r^2 0.25-0.4

###############################
## missingness
df <- read.table("both_sex_additive_height.1.profile",header=TRUE)
df <- read.table("EUR_male_height_prs_2.1.profile",header=TRUE)
df <- read.table("missing_height.vmiss",header=FALSE)
df <- df[df$CNT !=0,]
head(df)
x <- 1:nrow(df)
plot(x,df$CNT,pch=".",ylab="CNT",main="both-sex additive - height - p<1")
hist(df$CNT,xlab="CNT",ylab="frequency",main="both-sex additive - height - p<1")
hist(df$V5,xlab="variant missingness rate",ylab = "frequency", 
     main="height test set - variant missingness")
max(df$CNT) - min(df$CNT)


##############
df <- read.table("height_pheno_covar.txt",header=TRUE)
df_f <- df[df$sex==0,]; df_m = df[df$sex ==1,]
plot(df_f$height)
plot(df_m$height)


#######
sample <- sample(c(TRUE, FALSE), nrow(mtcars), replace=TRUE, prob=c(0.7,0.3))
train <- mtcars[sample, ]
test <- mtcars[!sample, ]
head(test)
