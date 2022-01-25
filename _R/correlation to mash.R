library(ggplot2)
library(ggrepel)
library(ggsci)
library(reshape2)
library(ggbreak)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)

# get correlation
setwd("~/Research/GWAS-frontera/LDSC")
ldsc <- read.csv("ldsc_results.txt", sep="\t")
ldsc <- ldsc[ldsc$Sex == "both_sex", c(1,2,6)]

# get mash weights
setwd("~/Research/GWAS-frontera/OLD/GWAS_Results_OLD/")
mash <- read.csv("sum_mash_weights.txt", sep="\t")

df <- merge(ldsc, mash, by.x = "Code", by.y = "phenotype")
df <- df %>%
  mutate(diff_weight = sum_weight_m - sum_weight_f) %>%
  mutate(nontrivial = 1 - (trivial/100)) %>%
  #select(c(2,3,8,9)) %>%
  mutate(facet = ifelse(Phenotype == "Testosterone", 1, 2))


# x axis is M>F weight amplification
# y axis is (1 - trivial weight)  trivial weight = perfect correlation and equal magnitude
ggplot(data= df, aes(x= Correlation, y= nontrivial)) +
  geom_point(size=3) +
  facet_grid(~facet, scales = "free_x", space = "free") +
  theme(strip.text.x = element_blank()) +
  labs(title = "Non-Trivial Weight on Covariance Matrice \nCompared with Genetic Correlation", x = "Genetic Correlation", 
       y = "Proportion of Non-trivial Weight \non Hypothesis Covariance Matrices") +
  geom_text_repel(aes(label = Phenotype), box.padding = 0.7, max.overlaps = Inf, size=3.5)


ggplot(data= df, aes(x= diff_weight, y= nontrivial)) +
  geom_point(size=3) +
  theme(strip.text.x = element_blank()) +
  labs(title = "Nontrivial Covariance Matrice Weight by Amplfication Signal", x = "Difference between fractions of variants with M>F and M<F effect", 
       y = "Proportion of Non-trivial Weight \non Hypothesis Covariance Matrices") +
  geom_text_repel(aes(label = Phenotype), box.padding = 0.7, max.overlaps = Inf, size=3.5)


