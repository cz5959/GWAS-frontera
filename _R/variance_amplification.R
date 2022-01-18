
# get nice phenotype names
setwd("~/Research/GWAS-frontera/LDSC/")
ldsc_df <- read.csv("ldsc_results.txt", sep="\t", colClasses = c(rep("character",2), rep("NULL",5)))
ldsc_df <- unique(ldsc_df)

# get phenotypic variance by sex
setwd("~/Research/GWAS-frontera/Phenotypes")
pheno_vars <- read.csv("pheno_var.txt", sep="\t")

# get mash weights
setwd("~/Research/GWAS-frontera/OLD/GWAS_Results_OLD/")
mash_weights <- read.csv("sum_mash_weights.txt", sep="\t")

# create column of ratio of male to female phenotypic variance
pheno_vars <- pheno_vars %>%
  merge(ldsc_df, by.x="pheno", by.y="Code") %>%
  mutate(ratio_mf = m_var / f_var) %>%
  select(c(1,6,7))

# create column of difference of male and female mash weights
mash_weights <- mash_weights %>%
  mutate(diff_mf = sum_weight_m - sum_weight_f) %>%
  select(c(1,4))

# combine phenotypic variance and mash weights df
df <- merge(pheno_vars, mash_weights, by.x = "pheno", by.y = "phenotype")

# split outliers (arm_fatfree_mass and testosterone)
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()

# scatter plot
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,40,20), limits=c(-70,50), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) + 
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x.bottom = element_blank())

g2 <- ggplot(data=df3, aes(x= diff_mf, y= ratio_mf)) +
  geom_point(size=2, color='#1d47a1') + 
  geom_hline(yintercept = 27.6) + geom_vline(xintercept = 95) +
  scale_x_continuous(breaks=c(85,90), limits=c(80,95), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text_repel(aes(label=Phenotype))

g3 <- ggplot(data=df1, aes(x= diff_mf, y= ratio_mf)) +
  geom_point(size=2, color='#1d47a1') +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,40,20), limits=c(-70,50), expand=c(0,0)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none") +
  geom_text_repel(aes(label=Phenotype), max.overlaps = Inf, box.padding=0.7)

g4 <- ggplot(df2, aes(x= diff_mf, y= ratio_mf)) + 
  geom_point(size=2, color='#1d47a1') +
  geom_vline(xintercept = 95) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=c(85,90), limits=c(80,95), expand=c(0,0)) +
  scale_y_continuous(breaks=c(1,2,3), limits=c(0,3.5), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) +
  geom_text_repel(aes(label=Phenotype), box.padding=1)

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g3)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, gB, g4, ncol=2, nrow=2, widths=c(5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Difference between fractions of variants with M>F and M<F effect", size=14),
                left = text_grob("Ratio of Male to Female Trait Variance", size=14, rot=90),
                top = text_grob("Amplification strongly correlates with phenotypic variance", face = "bold", size = 14))


