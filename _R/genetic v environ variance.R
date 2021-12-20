library(ggplot2)
library(ggrepel)
library(ggsci)
library(reshape2)
library(ggbreak)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)
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
  m <- df_pheno[df_pheno$sex == 1, 2] ; f <- df_pheno[df_pheno$sex == 0, 2]
  male_var <- var(m) ; male_varse <- sqrt((2*male_var^2) / (length(m)-1))
  female_var <- var(f) ; female_varse <- sqrt((2*female_var^2) / (length(f)-1))
  pheno_var <- rbind(pheno_var, data.frame(pheno=pheno, m_var=male_var, f_var=female_var,
                                           m_varse=male_varse, f_varse=female_varse))
}

# get heritability
setwd("~/Research/GWAS-frontera/LDSC")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t")
ldsc_df <- ldsc_df[c(1,3,4,5)] ; ldsc_df <- ldsc_df[ldsc_df$Sex != 'both_sex',]
ldsc_se <- dcast(ldsc_df, Code ~ Sex, value.var = c("h2.std.error"))
ldsc_df <- dcast(ldsc_df, Code ~ Sex, value.var = c("Heritability"))
colnames(ldsc_df) <- c("pheno","f_h2","m_h2") ; colnames(ldsc_se) <- c("pheno","f_h2_se","m_h2_se")
ldsc_df <- ldsc_df %>%
  mutate(f_h2e = 1-f_h2) %>%
  mutate(m_h2e = 1-m_h2)

# merge variance and heritability dataframe
df <- merge(merge(pheno_var, ldsc_df, by='pheno'), ldsc_se, by='pheno')
df <- df %>%
  mutate(geno_var_m = m_h2*m_var) %>% mutate(geno_var_f = f_h2*f_var) %>%
  mutate(env_var_m = m_h2e*m_var) %>% mutate(env_var_f = f_h2e*f_var) %>%
  mutate(geno_se_m = sqrt((m_varse^2*m_h2_se^2)+(m_h2^2*m_varse^2)+(m_var^2*m_h2_se^2))) %>% 
  mutate(geno_se_f = sqrt((f_varse^2*f_h2_se^2)+(f_h2^2*f_varse^2)+(f_var^2*f_h2_se^2))) %>%
  mutate(env_se_m = sqrt((m_varse^2*m_h2_se^2)+(m_h2e^2*m_varse^2)+(m_var^2*m_h2_se^2))) %>% 
  mutate(env_se_f = sqrt((f_varse^2*f_h2_se^2)+(f_h2e^2*f_varse^2)+(f_var^2*f_h2_se^2))) %>%
  # calculate ratio
  mutate(geno_var_ratio = geno_var_m/geno_var_f) %>% 
  mutate(env_var_ratio = env_var_m/env_var_f) %>%
  # estimate standard error of ratio
  mutate(geno_var_ratio_se = 
           (geno_var_ratio)*sqrt((geno_se_m^2/geno_var_m^2)+(geno_se_f^2/geno_var_f^2)) ) %>%
  mutate(env_var_ratio_se = 
           (env_var_ratio)*sqrt((env_se_m^2/env_var_m^2)+(env_se_f^2/env_var_f^2))) %>%
  select(c(1,20,21,22,23))

# ggplot
df$facet_y[df$env_var_ratio <= 10] <- 2
df$facet_y[df$env_var_ratio > 10] <- 1
df$facet_x[df$geno_var_ratio <= 2] <- 1
df$facet_x[df$geno_var_ratio > 3] <- 2
df$facet_x[df$geno_var_ratio > 10] <- 3

# plot ### OLD
# g <- ggplot(data=df, aes(x= geno_var_ratio, y= env_var_ratio)) +
#   geom_point(size=2, alpha=0.5) +
#   geom_abline(slope = 1, intercept = 0) + 
#   geom_errorbar(aes(ymin=env_var_ratio-env_var_ratio_se, ymax=env_var_ratio+env_var_ratio_se), 
#                 width=0, position="dodge") +
#   geom_errorbarh(aes(xmin=geno_var_ratio-env_var_ratio_se, xmax=geno_var_ratio+env_var_ratio_se), 
#                 position="dodge") +
#   labs(title="Genetic Variance Ratio to \nEnvironmental Variance Ratio", 
#   x="Genetic Variance Ratio", y="Environmental Variance Ratio") +
#   theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=15)) +
#   facet_grid(facet_y~facet_x, scales="free", space="free") +
#   geom_text_repel(aes(label=pheno), size=3, max.overlaps=100, segment.alpha=0.5) +
#   theme_pubclean() + scale_color_npg()
# g



## GRID EXTRA PLOTTING
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df1 <- mutate(df1, on = ifelse(abs(env_var_ratio-geno_var_ratio) <= env_var_ratio_se | 
                                 abs(geno_var_ratio-env_var_ratio) <= geno_var_ratio_se, 2,1))
df1_1 <- df1[(df1$env_var_ratio/df1$geno_var_ratio) >=1,]
df1_2 <- df1[df1$env_var_ratio/df1$geno_var_ratio <1,]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()

g1 <- ggplot(empty_df) + geom_point() + 
  scale_y_continuous(breaks=c(24.116), limits=c(23.45,24.78)) + 
  scale_x_continuous(breaks=c(0.5,1,1.5), limits=c(0.25,1.75)) +
  theme(axis.text.x=element_blank(), axis.title=element_blank(), axis.ticks.x=element_blank())

g2 <- ggplot(empty_df) + geom_point() +
  scale_x_continuous(breaks=c(3,3.4,3.8), limits=c(3,3.9)) +
  scale_y_continuous(breaks=c(24.116), limits=c(23.45,24.78)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank() )

g3 <- ggplot(data=df3, aes(x= geno_var_ratio, y= env_var_ratio)) +
  geom_point(size=2, alpha=0.5) + 
  geom_errorbar(aes(y=env_var_ratio, ymin=env_var_ratio-env_var_ratio_se, ymax=env_var_ratio+env_var_ratio_se), width=0) +
  geom_errorbarh(aes(y=env_var_ratio, xmin=geno_var_ratio-geno_var_ratio_se, xmax=geno_var_ratio+geno_var_ratio_se), height=0) +
  scale_y_continuous(breaks=c(24.116), limits=c(23.4,24.8)) +
  scale_x_continuous(breaks=c(77.715), limits=c(64.35,91.07)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  geom_text_repel(aes(label=pheno), size=3, segment.alpha=0.5)
  
g4 <- ggplot(data=df1, aes(x= geno_var_ratio, y= env_var_ratio, color=on)) +
  geom_point(size=2, alpha=0.5) + geom_abline(slope = 1, intercept = 0, color='red') +
  geom_errorbar(aes(y=env_var_ratio, ymin=env_var_ratio-env_var_ratio_se, ymax=env_var_ratio+env_var_ratio_se), 
                width=0, position="dodge") +
  geom_errorbarh(aes(y=env_var_ratio, xmin=geno_var_ratio-geno_var_ratio_se, xmax=geno_var_ratio+geno_var_ratio_se)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5)) +
  scale_x_continuous(breaks=c(0.5,1,1.5), limits=c(0.25,1.75)) +
  theme(axis.text.x=element_text(), axis.title=element_blank(), legend.position = "none") +
  geom_text_repel(data=(df1_1), aes(label=pheno), size=3, max.overlaps=Inf, segment.alpha=0.5, 
                  ylim=c(2.2,NA)) +
  geom_text_repel(data=(df1_2), aes(label=pheno), size=3, max.overlaps=Inf, segment.alpha=0.5, 
                  ylim=c(NA,0.4))

g5 <- ggplot(data=df2, aes(x= geno_var_ratio, y= env_var_ratio)) +
  geom_point(size=2, alpha=0.5) + geom_abline(slope = 1, intercept = 0, color='red') +
  geom_errorbar(aes(y=env_var_ratio, ymin=env_var_ratio-env_var_ratio_se, 
                    ymax=env_var_ratio+env_var_ratio_se), width=0) +
  geom_errorbarh(aes(y=env_var_ratio, xmin=geno_var_ratio-geno_var_ratio_se, xmax=geno_var_ratio+geno_var_ratio_se), height=0) +
  scale_x_continuous(breaks=c(3,3.4,3.8), limits=c(3,3.9)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5)) +
  theme(axis.text.y=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank()) +
  geom_text_repel(aes(label=pheno), size=3, max.overlaps=Inf, segment.alpha=0.5,
                  ylim=c(0,2.5))

g6 <- ggplot(empty_df) + geom_point() +
  scale_x_continuous(breaks=c(77.715), limits=c(77.70,77.73)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5)) +
  theme(axis.text.y=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank())

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g4)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, g3, gB, g5, g6, ncol=3, nrow=2, widths=c(6,2,1), heights=c(1,4))

annotate_figure(plot, 
                bottom = text_grob("Genetic Variance Ratio", size=10),
                left = text_grob("Environmental Variance Ratio", size=10, rot=90),
                top = text_grob("Genetic Variance Ratio to \n Environmental Variance Ratio", face = "bold", size = 14))

head(mtcars[,3])
mtcars$disp[3]
a <- "test"
mtcars[a] <- 0 
