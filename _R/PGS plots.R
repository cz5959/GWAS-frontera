
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)

# load files
setwd("~/Research/GWAS-frontera/GWAS_Results")
df_pgs <- read.csv("pgs_linear_results.txt",sep="\t")
setwd("~/Research/GWAS-frontera/LDSC")
df_ldsc <- read.csv("ldsc_results.txt",sep="\t")

# constrain phenotypes
df_pgs <- df_pgs[df_pgs$Phenotype %in% df_ldsc$Phenotype,]
df_ldsc <- df_ldsc[df_ldsc$Phenotype %in% df_pgs$Phenotype,]

# merge dataframes
df_pgs$Sex[df_pgs$Sex == 'combined'] <- 'both_sex'
df_scatter <- merge(df_pgs,df_ldsc,by=c("Phenotype","Sex"))

# model ratios
df_scatter$addboth_by_mash <- df_scatter$addboth_inc_r2 / df_scatter$mash_inc_r2
df_scatter$addboth_by_add <- df_scatter$addboth_inc_r2 / df_scatter$add_inc_r2
df_scatter$mash_by_add <- df_scatter$mash_inc_r2 / df_scatter$add_inc_r2

# both_sex, female and male dataframe
df_scatter_bs <- df_scatter[df_scatter$Sex == 'both_sex',c(1,2,5,9,13,15,17)]
df_scatter_f <- df_scatter[df_scatter$Sex == 'female',c(1,2,5,9,13,15,17)]
df_scatter_m <- df_scatter[df_scatter$Sex == 'male',c(1,2,5,9,13,15,17)]
df_scatter_2 <- df_scatter[c(1,2,6,10,14,15,17)]

# melt
df_scatter_long <- melt(df_scatter_2, id.vars=c("Phenotype","Sex","Heritability","Correlation"))
df_scatter_bs_long <- melt(df_scatter_bs, id.vars=c("Phenotype","Sex","Heritability","Correlation"))
df_scatter_f_long <- melt(df_scatter_f, id.vars=c("Phenotype","Sex","Heritability","Correlation"))
df_scatter_m_long <- melt(df_scatter_m, id.vars=c("Phenotype","Sex","Heritability","Correlation"))

# heritability ratios
h2f_by_h2m <- df_scatter_f$Heritability / df_scatter_m$Heritability
h2b_by_h2m <- df_scatter_bs$Heritability / df_scatter_m$Heritability
h2b_by_h2f <- df_scatter_bs$Heritability / df_scatter_f$Heritability
df_scatter$h2f_by_h2m <- rep(h2f_by_h2m, each = 3)
df_scatter$h2b_by_h2m <- rep(h2b_by_h2m, each = 3)
df_scatter$h2b_by_h2f <- rep(h2b_by_h2f, each = 3)

df_scatter_3 <- df_scatter[c(1,2,19,20,21,22,23,24)]
df_scatter_long_3 <- melt(df_scatter_3, id.vars=c("Phenotype","Sex","h2f_by_h2m","h2b_by_h2m", "h2b_by_h2f"))
head(df_scatter_long_3)

### PLOTS ###

# r2 by heritability
ggplot(df_scatter, aes(x=Heritability , y=mash_inc_r2)) +
  geom_point(size=2) +
  geom_smooth(method="lm", alpha=0.2) +
  ylab("mash Incremental R2") +
  ggtitle("mash Incremental R2 by Heritability") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=12), legend.text=element_text(size=12)) +
  facet_grid(~Sex)

# r2 by heritability by sex
ggplot(df_scatter_bs_long, aes(x=Heritability , y=value, color=variable)) +
  geom_point(size=2) +
  geom_smooth(method="lm", alpha=0.2) + 
  labs(title="R2 by Heritability - Both-Sex", x="Heritability", y="R2", color="R2") +
  scale_color_hue(labels = c("additive both-sex", "additive same-sex", "mash")) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=12), legend.text=element_text(size=12))

# r2 by heritability
df_scatter_long_2 <- df_scatter_long[df_scatter_long$Phenotype != 'Testosterone',]
labels <- c(addboth_inc_r2="Additive Both-Sex", add_inc_r2="Additive Same-Sex", mash_inc_r2="mash")
ggplot(df_scatter_long, aes(x=Heritability , y=value, color=variable)) +
  geom_point(size=2) +
  geom_smooth(method="lm", alpha=0.2) + 
  labs(title="Incremental R2 by Heritability", x="Heritability", y="Incremental R2") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=12), legend.text=element_text(size=12)) +
  facet_grid(~Sex, labeller=labeller(variable=labels))

# add both r2 / mash r2  by heritability
df_scatter_omit <- df_scatter[!df_scatter$Phenotype %in% c('Testosterone'),]

ggplot(df_scatter_omit, aes(x=h2f_by_h2m, y=addboth_by_mash)) +
  geom_point(size=2) +
  geom_smooth(method="lm",alpha=0.1) +
  labs(title="Ratio of Additive Both-Sex inc. R2 to mash inc. R2", x="Ratio of Female h2 to Male h2", y="Ratio of Additive Both-Sex to mash inc. R2") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=12), legend.text=element_text(size=12)) +
  facet_grid(~Sex) + 
  geom_text_repel(label=df_scatter_omit$Phenotype, size=3, max.overlaps=100, segment.alpha=0.5)

# ratio to ratio
ggplot(df_scatter_long_3, aes(x=h2f_by_h2m, y=value, color=variable)) +
  geom_point(size=2) +
  geom_smooth(method="lm", alpha=0.15, aes(fill=variable), show.legend=FALSE) +
  labs(title="Comparing Model Ratios against Heritability Ratio", x="Ratio of Female h2 to Male h2", y="Model to Model Ratio", col="Model Ratios") +
  scale_color_hue(labels = c("add. both-sex / mash", "add. both-sex / add. same-sex", "mash / add. same-sex")) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=12), legend.text=element_text(size=12)) +
  facet_grid(~Sex) +
  stat_cor(method='pearson', p.accuracy=0.001)

