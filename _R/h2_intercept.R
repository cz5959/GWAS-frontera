require("dplyr")
require("tidyr")
require("ggplot2")
require("grid")
require("gridExtra")
require("ggpubr")

setwd("~/Research/Phenotypes/")
names <- read.csv("pheno_names.txt", sep="\t")
setwd("~/Research/GWAS-frontera/mash")
mash <- read.csv("mash_weights.txt", sep="\t")
setwd("~/Research/GWAS-frontera/LDSC")
df <- read.csv("h2_intercept.txt", sep="\t")

df <- merge(merge(df, names, by="Phenotype"), mash, by.x="Code", by.y="phenotype")

df <- df %>%
  filter(Sex != "both_sex") %>% 
  mutate(f_m_diff = sum_weight_f - sum_weight_m)
nrow(df)  

rects <- data.frame(ystart = seq(0.5,26.5,1), yend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
ggplot(df, aes(x=intercept, y=reorder(Phenotype, f_m_diff), col=Sex)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(y=Phenotype, xmin=intercept-intercept.std.error, xmax=intercept+intercept.std.error), 
                 height=0, position=position_dodge(width=0.5)) +
  geom_rect(data=rects, aes(ymin=ystart,ymax=yend,xmin=1,xmax=1.3), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text = element_text(size=10), 
        axis.title.y = element_blank(), axis.title.x = element_text(size=12),
        plot.margin = margin(5,10,5.5,5.5)) +
  scale_color_manual(labels = c("female", "male"), values=c("#d67629","#207335")) +
  scale_x_continuous(breaks=seq(1,1.3,0.05),expand=c(0,0))
  
df[df$Phenotype == "Weight",]
