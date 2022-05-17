library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)
library(ggrepel)

pheno_list <- c("height","bmi","RBC_count","IGF1","creatinine","weight","calcium",
                "protein_total","whole_body_fat_mass","urate","arm_fatfree_mass_L",
                "arm_fatfree_mass_R", "eosinophil_perc", "lymphocyte_perc", "waist_circ",
                "hip_circ", "waist_to_hip", "wth_bmi_adj","diastolicBP_auto","systolicBP_auto",
                "albumin", "pulse_rate", "urea", "SHBG", "FVC_best", "HbA1c")

# get lm results from each testosterone bin
bin_fun <- function(data, n, sex) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  for (i in 1:n) {
    # linear regression
    bin <- data[cuts == levels(cuts)[i],]
    model <- lm(paste0("pheno ~ SCORE"), data = bin)
    beta <- model$coefficients[2]
    T_mean <- mean(bin$testosterone)
    age_mean <- mean(bin$age)
    results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Age = age_mean, Sex=sex))
  }
  return(results)
}

get_corr <- function(data) {
  # correlation between Beta and testosterone
  corrs_row <- NULL
  for (sex in c("male", "female")) {
    data_sub <- data[data$Sex == sex,]
    model <- lm(Beta ~ Age, data_sub)
    corr <- cor.test(data_sub$Testosterone, residuals(model))
    corr_est <- corr$estimate
    corr_err <- corr$conf.int[2] - corr_est
    corrs_row <- rbind(corrs_row, data.frame(Pheno=pheno, Est=corr_est, Err=corr_err, Sex=sex))
  }
  corrs_row$est_diff <- abs(corrs_row[1,2] - corrs_row[2,2])
  corrs_row$err_sum <- abs(corrs_row[1,3] + corrs_row[2,3])
  return(corrs_row)
}


corrs_result <- NULL
for (pheno in pheno_list) {
  print(pheno)
  # phenotype
  setwd("~/Research/Phenotypes")
  df_testosterone <- read.csv("pheno_testosterone.txt", sep="\t", colClasses = c("NULL","integer","numeric"))
  df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
  df_sex <- read.csv("sex_ids.txt", sep="\t")
  df_age <- read.csv("covariates.txt", sep="\t", colClasses =c("NULL", "integer", rep("NULL",11), "integer", rep("NULL", 3)))
  colnames(df_age) <- c('IID','age')
  colnames(df_pheno) <- c('IID','pheno')
  
  # PGS scores
  setwd(paste0("~/Research/GWAS-frontera/GWAS_results/",pheno))
  file_name <- list.files(pattern="both_sex_additive_")
  df_both <- read.csv(file_name,sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))

  # merge dataframes
  df <- merge(merge(merge(df_testosterone, df_sex, by='IID'), df_age, by = 'IID'), df_pheno, by='IID')
  df <- merge(df,df_both,by='IID')

  # order by testosterone
  df <- df[order(df$testosterone),]
  # label then split by sex
  df$sex[df$sex == 1] <- 'male'
  df$sex[df$sex == 0] <- 'female'
  df_m <- df[df$sex == 'male',]
  df_f <- df[df$sex == 'female',]

  # call function to get betas for each bin
  m_results <- bin_fun(df_m,10,'male')
  f_results <- bin_fun(df_f,10,'female')
  results <- rbind(m_results, f_results)

  # correlation between Beta and testosterone
  corrs_result <- rbind(corrs_result, get_corr(results))
}
#setwd("~/Research/Phenotypes")
#write.table(corrs_result, file="G_corr_testosterone_age.txt", sep="\t", row.names=FALSE)
corrs_result <- read.csv("G_corr_testosterone_age.txt", sep="\t")
# get nice phenotype names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)
corrs_result <- merge(corrs_result, ldsc_df, by.x="Pheno", by.y="Code")

## plot
rects <- data.frame(xstart = seq(0.5,26.5,1), xend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
rects <- rects[1:26,]
ggplot(corrs_result, aes(x=reorder(Phenotype, est_diff), y=Est, color=Sex)) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=Est-(1.645*Err), ymax=Est+(1.645*Err)), 
                alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  geom_rect(data=rects, aes(xmin=xstart,xmax=xend,ymin=-1.5,ymax=1.2), 
            inherit.aes = FALSE, alpha=0.2, fill = c(rep(c("grey","white"),13))) +
  labs(x="", y="R") +
  scale_y_continuous(expand=c(0.01,0)) +
  theme_classic() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_blank(),
        legend.position = "top", legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  scale_color_manual(values=c("#d67629","#207335")) +
  coord_flip()

head(corrs_result)
