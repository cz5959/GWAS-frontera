
# get nice phenotype names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)

# get phenotypic variance by sex
setwd("~/Research/GWAS-frontera/Phenotypes")
pheno_vars <- read.csv("pheno_var.txt", sep="\t")

# get mash weights
setwd("~/Research/GWAS-frontera/mash/")
mash_weights <- read.csv("mash_weights.txt", sep="\t")

# create column of ratio of male to female phenotypic variance
pheno_vars <- pheno_vars %>%
  merge(ldsc_df, by.x="pheno", by.y="Code") %>%
  mutate(ratio_mf = m_var / f_var) %>%
  select(c(1,6,7))

# create column of difference of male and female mash weights
mash_weights <- mash_weights %>%
  mutate(diff_mf = sum_weight_m - sum_weight_f) %>%
  select(c(1,6))

# combine phenotypic variance and mash weights df
df <- merge(pheno_vars, mash_weights, by.x = "pheno", by.y = "phenotype")

df_corr <- df[df$pheno != 'testosterone',]
model <- lm(ratio_mf ~ diff_mf, df_corr)
model <- cor.test(df_corr$diff_mf, df_corr$ratio_mf)
summary(model)
model

# split outliers (arm_fatfree_mass and testosterone)
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()
anno_size=2.8
# scatter plot
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) + 
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=10),
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
  geom_point(size=2, color='#2b62d9') +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  geom_abline(slope = 0.014007, intercept = 1.220738, color="gray", size=0.8) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=10)) +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)

g4 <- ggplot(df2, aes(x= diff_mf, y= ratio_mf)) + 
  geom_point(size=2, color='#2b62d9') +
  geom_vline(xintercept = 100) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_abline(slope = 0.014007, intercept = 1.220738, color="gray", size=0.8) +
  scale_x_continuous(breaks=c(85,90,95), limits=c(80,100), expand=c(0,0)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(size=10),
        axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size) + 
  coord_cartesian(clip="off")

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g3)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, gB, g4, ncol=2, nrow=2, widths=c(5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Difference between fractions of variants with M>F and M<F effect", size=12),
                left = text_grob("Ratio of Male to Female Phenotypic Variance", size=12, rot=90),
                top = text_grob("Phenotypic Variance Strongly Correlated with Amplification", size = 16))


