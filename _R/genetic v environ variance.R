library(ggplot2)
library(ggrepel)
library(ggsci)
library(reshape2)
library(ggbreak)

# get phenotypic variance
setwd("~/Research/GWAS-frontera/Phenotypes")
pheno_var <- NULL
pheno_list <- c("height","bmi","testosterone","RBC_count","IGF1","creatinine","weight","calcium",
                "protein_total","whole_body_fat_mass","urate","arm_fatfree_mass_L",
                "arm_fatfree_mass_R", "eosinophil_perc", "lymphocyte_perc", "waist_circ",
                "hip_circ", "waist_to_hip", "wth_bmi_adj","diastolicBP_auto","systolicBP_auto",
                "albumin", "pulse_rate")
df_sex <- read.csv("sex_ids.txt", sep="\t")
for (pheno in pheno_list) {
  df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
  df_pheno <- merge(df_pheno, df_sex, by='IID')
  male_var <- var(df_pheno[df_pheno$sex == 1, 2])
  female_var <- var(df_pheno[df_pheno$sex == 0, 2])
  pheno_var <- rbind(pheno_var, data.frame(pheno=pheno, m_var=male_var, f_var=female_var))
}

# get heritability
setwd("~/Research/GWAS-frontera/LDSC")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t")
ldsc_df <- ldsc_df[c(1,3,4)]
ldsc_df <- ldsc_df[ldsc_df$Sex != 'both_sex',]
ldsc_df <- dcast(ldsc_df, Code ~ Sex, value.var = "Heritability")
colnames(ldsc_df) <- c("pheno","f_h2","m_h2")
ldsc_df$f_h2e <- 1 - ldsc_df$f_h2
ldsc_df$m_h2e <- 1 - ldsc_df$m_h2

# merge variance and heritability dataframe
df <- merge(pheno_var, ldsc_df, by='pheno')
df$geno_var_m <- df$m_h2 * df$m_var ; df$geno_var_f <- df$f_h2 * df$f_var
df$env_var_m <- df$m_h2e * df$m_var ; df$env_var_f <- df$f_h2e * df$f_var
df$geno_var_ratio <- df$geno_var_m / df$geno_var_f
df$env_var_ratio <- df$env_var_m / df$env_var_f
df <- df[c(1,12,13)]

df$facet_y[df$env_var_ratio <= 2.5] <- 1
df$facet_y[df$env_var_ratio > 2.5] <- 2
df$facet_y[df$env_var_ratio > 20] <- 3
df$facet_x[df$geno_var_ratio <= 10] <- 2
df$facet_x[df$geno_var_ratio > 10] <- 1

# plot
g <- ggplot(data=df, aes(x= geno_var_ratio, y= env_var_ratio)) +
  geom_point(size=2) +
  labs(title="Genetic Variance Ratio to \nEnvironmental Variance Ratio", 
  x="Genetic Variance Ratio", y="Environmental Variance Ratio") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=15)) +
  facet_grid(facet_x~facet_y, scales="free", space="free") +
  geom_text_repel(label=df$pheno, size=3, max.overlaps=100, segment.alpha=0.5) +
  theme_bw() + scale_color_npg()
g




