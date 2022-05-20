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

# get nice phenotype names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)

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
    stderror <- summary(model)$coefficients[2,2]
    T_mean <- mean(bin$testosterone)
    results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Error=stderror, Sex=sex))
  }
  return(results)
}

get_corr <- function(data) {
  # correlation between Beta and testosterone
  corrs_row <- NULL
  for (sex in c("male", "female")) {
    data_sub <- data[data$Sex == sex,]
    corr <- cor.test(data_sub$Testosterone, data_sub$Beta)
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
  setwd("~/Research/GWAS-frontera/Phenotypes")
  df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
  df_sex <- read.csv("sex_ids.txt", sep="\t")

  # PGS scores
  setwd(paste0("~/Research/GWAS-frontera/GWAS_results/",pheno))
  file_name <- list.files(pattern="both_sex_additive_")
  df_both <- read.csv(file_name,
                      sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  # testosterone pgs scores
  setwd("~/Research/GWAS-frontera/GWAS_results/testosterone")
  #df_testosterone <- read.csv("both_sex_additive_testosterone.1e-5.profile",
  #                            sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  ## sex specific testosterone
  df_T_f <- read.csv("female_additive_testosterone.1e-8.profile", 
                     sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  df_T_m <- read.csv("male_additive_testosterone.1e-5.profile", 
                     sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  df_testosterone <- merge(df_T_f,df_T_m, by="IID")

  # merge dataframes
  df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
  ## sex specific testosterone
  df$SCORE.x <- ifelse(df$sex == 1, df$SCORE.x, df$SCORE.y)   
  df <- df[-c(3)]
  
  df <- merge(df,df_both,by='IID')
  colnames(df) <- c("IID", "testosterone", "sex", "pheno", "SCORE")

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

# save corr and slope files
setwd("~/Research/Phenotypes")
#write.table(corrs_result, file="G_corr_testosterone_sexspecificpgs.txt", sep="\t", row.names=FALSE)
#corrs_result <- read.csv("G_corr_testosterone_pgs.txt", sep="\t")
corrs_result <- read.csv("G_corr_testosterone_sexspecificpgs.txt", sep="\t")

# ldsc correlation
setwd("~/Research/GWAS-frontera/LDSC")
ldsc_df <- read.csv("ldsc_results.txt",sep="\t")
ldsc_df <- ldsc_df[ldsc_df$Sex == 'both_sex',c(1,2,6)]

format_result <- function(result) {
  result <- result[order(result$Pheno),]
  result <- merge(result, ldsc_df, by.x='Pheno', by.y='Code')
  return(result)
}
corrs_result <- format_result(corrs_result)
head(corrs_result)
#####################################################
# new_pheno <- c("Albumin", "Arm fat free\nmass (L)", "Arm fat free\nmass (R)", "BMI", "Calcium", "Creatinine",
#   "Diastolic BP", "Eosinophil\npercentage", "Forced vital\ncapacity", "HbA1c", 
#   "Height", "Hip\ncircumference", "IGF-1", "Lymphocyte\npercentage", "Total protein", "Pulse rate", 
#   "RBC count", "SHBG", "Systolic BP", "Urate", "Urea", "Waist\ncircumference", "Waist to hip\nratio",
#   "Weight", "Whole body\nfat mass" , "Waist:hip\n(bmi adj.)") 
# corrs_result$Phenotype <- rep(new_pheno, each=2)

# corr plot scatter
# 90% confidence interval --> mean +/- 1.645*SE ; z-score=1.645
rects <- data.frame(xstart = seq(0.5,26.5,1), xend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
rects <- rects[1:26,]
plot <- ggplot(corrs_result, aes(x= reorder(Phenotype, est_diff), y=Est, color=Sex)) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=Est-(1.645*Err), ymax=Est+(1.645*Err)), alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  geom_rect(data=rects, aes(xmin=xstart,xmax=xend,ymin=-1.7,ymax=1.15), 
            inherit.aes = FALSE, alpha=0.2, fill = c(rep(c("grey","white"),13))) +
  labs(t="R - Correlation between Testosterone PGS and \nEffect of PGS on Phenotype") +
  scale_y_continuous(expand=c(0,0), breaks=seq(-1.5,1,0.5)) +
  theme_classic() + 
  theme(axis.text = element_text(size=10), axis.title = element_blank(), plot.title=element_blank(),
        legend.position = "none") +
  scale_color_manual(values=c("#d67629","#207335")) +
  coord_flip()
annotate_figure(plot, 
                bottom = text_grob("R - Correlation between Testosterone PGS and \nEffect of Polygenic Score on Phenotype", size=12))


