library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggrepel)

# get phenotypic mean and variance
setwd("~/Research/Phenotypes")
set.seed(1)
pheno_stats <- NULL
pheno_list <- c("height","bmi","testosterone","RBC_count","IGF1","creatinine","weight","calcium",
                "protein_total","whole_body_fat_mass","urate","arm_fatfree_mass_L",
                "arm_fatfree_mass_R", "eosinophil_perc", "lymphocyte_perc", "waist_circ",
                "hip_circ", "waist_to_hip", "wth_bmi_adj","diastolicBP_auto","systolicBP_auto",
                "albumin", "pulse_rate", "urea", "SHBG", "FVC_best", "HbA1c")

df_sex <- read.csv("sex_ids.txt", sep="\t")
for (pheno in pheno_list) {
   df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
   df_pheno <- merge(df_pheno, df_sex, by='IID')
   m <- df_pheno[df_pheno$sex == 1, 2] ; f <- df_pheno[df_pheno$sex == 0, 2]
   male_var <- var(m)
   female_var <- var(f)
   male_mean <- mean(m)
   female_mean <- mean(f)

   # bootstrap
   # male_varse <- replicate(100,var(sample(m,replace=T)))
   # female_varse <- replicate(100,var(sample(f,replace=T)))
   # male_varse <- sqrt( sum((male_varse - mean(male_varse))^2) / (100-1) )
   # female_varse <- sqrt( sum((female_varse - mean(female_varse))^2) / (100-1) )
   male_varse <- mean(replicate(100, sd(sample(m,replace=T))/sqrt(length(m))))
   female_varse <- mean(replicate(100, sd(sample(f,replace=T))/sqrt(length(f))))

   pheno_stats <- rbind(pheno_stats, data.frame(pheno=pheno, m_mean=male_mean, f_mean=female_mean, m_var=male_var, f_var=female_var,
                                            m_varse=male_varse, f_varse=female_varse))
}
head(pheno_stats)
write.table(pheno_stats, file="pheno_meanvar.txt", sep="\t", row.names=FALSE, quote=FALSE)

setwd("~/Research/Phenotypes")
pheno_var <- read.csv("pheno_meanvar.txt", sep="\t")
# get nice names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)
df <- merge(ldsc_df, pheno_var, by.x="Code", by.y="pheno")

# pheno mean ratio and variance ratio M:F
df <- df %>%
  mutate(mean_ratio = m_mean / f_mean) %>%
  mutate(var_ratio = m_var / f_var) %>%
  mutate(mean_diff = m_mean - f_mean)

# PLOT
head(df)
df1 <- df[!df$Code %in% c("testosterone","wth_bmi_adj"),]
df2 <- df[df$Code == "testosterone",]
df3 <- df[df$Code == "wth_bmi_adj",]
empty_df <- data.frame()
anno_size=2.8
axis_text=9

# scatter plot
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  scale_x_continuous(breaks=c(-1.2,-1.1), limits=c(-1.25,-1.1), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=axis_text),
        axis.ticks.x=element_blank(), axis.line.x.bottom = element_blank())

g2 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 1, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(trans="log10", breaks=seq(0.6,1.5,0.1), limits=c(0.55,1.7), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line= element_blank())

g3 <- ggplot(data=df2, aes(x= mean_ratio, y= var_ratio)) +
  geom_point(size=2, color='#2b62d9') + 
  geom_hline(yintercept = 27.6) + geom_vline(xintercept = 11) +
  scale_x_continuous(breaks=c(10.4,10.8), limits=c(10,11), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

# row2

g4 <- ggplot(data=df3, aes(x= mean_ratio, y= var_ratio)) + 
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_point(size=2, color='#2b62d9') + 
  scale_x_continuous(breaks=c(-1.2,-1.1), limits=c(-1.25,-1.0), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3,seq(0.5,3.5,0.25)), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text)) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

g5 <- ggplot(data=df1, aes(x= mean_ratio, y= var_ratio)) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_smooth(method="lm", color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  scale_x_continuous(trans="log10", breaks=seq(0.6,1.5,0.1), limits=c(0.55,1.7), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text),
        axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)
g5

g6 <- ggplot(empty_df) + geom_point() +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 11) +
  geom_smooth(method="lm", data=df1, aes(mean_ratio, var_ratio), inherit.aes = F, color="gray", size=0.8, se= F, fullrange=T) +
  scale_x_continuous(breaks=c(10.4,10.8), limits=c(10,11), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=axis_text),
        axis.ticks.y=element_blank(), axis.line.y = element_blank())

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g4)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, g3, gB, g5, g6, ncol=3, nrow=2, widths=c(1.5,5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Ratio of Male to Female Phenotypic Mean", size=11),
                left = text_grob("Ratio of Male to Female Phenotypic Variance", size=11, rot=90))


# mean_diff
df1 <- df1[df1$Code != "urate",]
ggplot(data=df1, aes(x= mean_diff, y= var_ratio)) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_smooth(method="lm", color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  theme_classic2() +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)
