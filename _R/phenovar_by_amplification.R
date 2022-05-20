
# get nice phenotype names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)

# get phenotypic variance by sex
setwd("~/Research/Phenotypes")
pheno_vars <- read.csv("pheno_meanvar.txt", sep="\t")

# get mash weights
setwd("~/Research/GWAS-frontera/mash/")
mash_weights <- read.csv("mash_weights.txt", sep="\t")

# create column of ratio of male to female phenotypic variance
pheno_vars <- pheno_vars %>%
  merge(ldsc_df, by.x="pheno", by.y="Code") %>%
  mutate(ratio_mf = m_var / f_var) %>%
  mutate(mean_ratio = m_mean / f_mean) %>%
  select(c(1,8,9,10))

# create column of difference of male and female mash weights
mash_weights <- mash_weights %>%
  mutate(diff_mf = sum_weight_m - sum_weight_f) %>%
  select(c(1,6))

# combine phenotypic variance and mash weights df
df <- merge(pheno_vars, mash_weights, by.x = "pheno", by.y = "phenotype")

df_corr <- df[df$pheno != 'testosterone',]
model <- cor.test(df_corr$diff_mf, df_corr$ratio_mf)
model
df_corr <- df_corr[df_corr$pheno != 'wth_bmi_adj',]
model <- cor.test(df_corr$mean_ratio, df_corr$diff_mf)
model
head(df_corr)
# split outliers (arm_fatfree_mass and testosterone)
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()
anno_size=2.8
axis_text = 9
# scatter plot 3.5x5.5
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=axis_text),
        axis.ticks.x=element_blank(), axis.line.x.bottom = element_blank())

g2 <- ggplot(data=df3, aes(x= diff_mf, y= ratio_mf)) +
  geom_point(size=2, color='#2b62d9') + 
  geom_hline(yintercept = 27.6) + geom_vline(xintercept = 100) +
  scale_x_continuous(breaks=c(85,90,95), limits=c(80,100), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

g3 <- ggplot(data=df1, aes(x= diff_mf, y= ratio_mf)) +
  geom_smooth(method="lm", data=df_corr, aes(x=diff_mf, y=ratio_mf), inherit.aes = F,
              color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3,seq(0.5,3.5,0.25)), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text)) +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)

g4 <- ggplot(df2, aes(x= diff_mf, y= ratio_mf)) + 
  geom_smooth(method="lm", data=df_corr, aes(x=diff_mf, y=ratio_mf), inherit.aes = F,
              color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  geom_vline(xintercept = 100) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=c(85,90,95), limits=c(80,100), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(size=axis_text),
        axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size) 

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g3)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, gB, g4, ncol=2, nrow=2, widths=c(5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Difference between fractions of variants with M>F and M<F effect", size=11),
                left = text_grob("Ratio of Male to Female Phenotypic Variances", size=12, rot=90),
                top = text_grob("Phenotypic Variance Strongly Correlated with Amplification", size = 14))


## 
setwd("~/Research/Phenotypes")
pheno_var <- read.csv("pheno_meanvar.txt", sep="\t")
pheno_var <- pheno_var %>%
  merge(ldsc_df, by.x="pheno", by.y="Code") %>%
  mutate(ratio_mf = m_mean / f_mean) %>%
  mutate(diff_mf = m_mean - f_mean)

df <- merge(pheno_var, mash_weights, by.x = "pheno", by.y = "phenotype")
df <- df[c(1,8,9,10,11)]
colnames(df) <- c("pheno", "Phenotype", "ratio_pheno","diff_pheno","diff_amplification")

ggplot(df1, aes(x=diff_pheno, y=diff_amplification)) +
  geom_point()

df1 <- df[!df$pheno %in% c("testosterone","wth_bmi_adj"),]
model <- cor.test(df1$ratio_pheno, df1$diff_amplification)
model
