
#######################################
pheno="asthma"
code=1111
library("crayon", lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")
library("dplyr", lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

# grab from non-cancer illness self-reported 
df <- read.table("NCI.txt",sep="\t",header=FALSE)
df_new <- df[, colSums(is.na(df[3:6])) < nrow(df)] # removed 38 columns 

# now just get phenotype based on code
ht <- df_new %>% filter_all(any_vars(. %in% code))

# 1 for controls; 2 for cases
df_new$pheno = rep(1,nrow(df_new))
df_new$pheno[df_new$V1 %in% ht$V1] <- 2
df <- df_new[c("V1","V2","pheno")]            ## TODO test
colnames(df) <- c("#FID","IID",pheno)

# save as txt file
write.table(df,file=paste0("pheno_",pheno,".txt"),sep="\t",row.names=FALSE,quote=FALSE)

#######################################
# change binary coding from 0/1 to 1/2
pheno="diabetes"
df <- read.table(paste0("pheno_",pheno,".txt"),sep="\t",header=FALSE)
df$V3[df$V3 == 1] <- 2
df$V3[df$V3 == 0] <- 1
df <- df[df$V3 > 0,]
colnames(df) <- c("#FID","IID",pheno)

write.table(df,file=paste0("pheno_",pheno,".txt"),sep="\t",row.names=FALSE,quote=FALSE)