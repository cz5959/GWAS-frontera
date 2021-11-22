library(ggplot2)

setwd("~/Research/GWAS-frontera/Phenotypes")
df_test <- read.csv("pheno_testosterone.txt", sep="\t")
df_cov <- read.csv("covariates.txt", sep="\t")

df_cov <- df_cov[c(2,13,14)]
df_test <- df_test[c(2,3)]

df <- merge(df_test,df_cov,by='IID')
colnames(df) <- c('IID', 'testosterone', 'sex', 'birthyear')
df$sex[df$sex == 1] <- 'male'
df$sex[df$sex == 0] <- 'female'

for (i in 1:nrow(df)) {
  if (df$birthyear < 1940) {
    df$group[i]
  }
}

ggplot(df, aes(x = testosterone, fill=sex)) +
  geom_histogram(bins=20, position="dodge")

