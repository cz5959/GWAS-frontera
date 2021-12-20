require("dplyr")
require("tidyr")
require("ggplot2")
require("gridExtra")
require("ggpubr")

setwd("~/Research/GWAS-frontera/LDSC")
df <- read.csv("ldsc_results.txt", sep="\t")

both_h2 <- rep(df[df$Sex == "both_sex", "Heritability"], each=3)

df <- df %>%
  mutate(relative_h2 = Heritability / both_h2) %>%
  mutate(relative_h2_se = ((Heritability + h2.std.error) / both_h2) - relative_h2 ) %>%
  arrange(Correlation) %>%
  mutate(Phenotype = factor(Phenotype, levels=unique(Phenotype)))

head(df)

ggplot(df, aes(x=relative_h2, y=Phenotype, col=factor(Sex))) +
  geom_point(size=2) +
  theme_pubclean() +
  theme(legend.position = "top", axis.text = element_text(size=11)) +
  scale_color_npg() +
  labs(title="Genetic Correlation by Relative Heritability", x="Relative Heritability", 
       y="Phenotype", col="Sex")
