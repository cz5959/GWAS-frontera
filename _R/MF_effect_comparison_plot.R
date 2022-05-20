library(ggplot2)
library(ggpubr)

setwd("~/Research/GWAS-frontera/mash/mash_posteriors")
df <- read.csv("mf_posteriors.txt", sep="\t")
head(df)

df <- df %>%
  mutate(slope_diff = posterior.slope..f.m. - raw.estimate.slope..f.m.)

rects <- data.frame(xstart = seq(0.5,27.5,1), xend = seq(1.5,28.5,1))
rects <- rects[1:27,]
ggplot(df, aes(x = reorder(phenotype,slope_diff), y = slope_diff, color=m.or.f.greater.var)) +
  geom_point() +
  labs(y="Difference between M:F posterior effect slope and \n M:F raw estimate slope") +
  theme_classic() +
  scale_color_manual(values=c("#d67629","#207335")) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=10), axis.text = element_text(size=8.5),
        legend.position = "none", plot.margin = margin(10,10,10,10)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_rect(data=rects, aes(xmin=xstart, xmax=xend, ymin=-0.01, ymax=0.5), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  geom_text(aes(label=alpha), hjust=-0.2, size=3) +
  coord_flip()

