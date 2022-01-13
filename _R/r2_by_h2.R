require("dplyr")
require("tidyr")
require("ggplot2")
require("grid")
require("gridExtra")
require("ggpubr")

setwd("~/Research/GWAS-frontera/LDSC")
df <- read.csv("ldsc_results.txt", sep="\t")

both_h2 <- rep(df[df$Sex == "both_sex", "Heritability"], each=3)

df <- df %>%
  mutate(relative_h2 = Heritability / both_h2) %>%
  mutate(relative_h2_se = ((Heritability + h2.std.error) / both_h2) - relative_h2 ) %>%
  mutate(Sex = factor(Sex, levels=c("female","both_sex","male"))) %>%
  arrange(Correlation, Phenotype, Sex) %>%
  mutate(Phenotype = factor(Phenotype, levels=unique(Phenotype))) 

rects <- data.frame(ystart = seq(0.5,26.5,1), yend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
p1 <- ggplot(df, aes(x=relative_h2, y=Phenotype, col=(Sex))) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(y=Phenotype, xmin=relative_h2-relative_h2_se, xmax=relative_h2+relative_h2_se), 
                 height=0, position=position_dodge(width=0.5)) +
  geom_rect(data=rects, aes(ymin=ystart,ymax=yend,xmin=0.5,xmax=2.5), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  scale_x_continuous(breaks=c(0.5,1,1.5,2), limits=c(0.5,2.5), expand=c(0,0)) +
  theme_bw() +
  #theme_classic() +
  theme(legend.position = "top", legend.text = element_text(size=11),
        axis.text = element_text(size=11), axis.text.y = element_text(hjust=0), 
        axis.title.y = element_blank(), axis.title.x = element_text(size=12),
        plot.margin = margin(5.5,5.5,5.5,0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(labels = c("female", "both sex", "male"), values=c("#d67629","#1d47a1","#207335")) +
  labs(x="Relative Heritability", y="Phenotype", col="Sex")

corr_df <- df %>%
  select(c(2,6)) %>%
  distinct() %>%
  mutate(pheno_point = seq(-0.03,1.03,by=1.06/(length(Phenotype)-1)))

p2 <- ggplot(corr_df, aes(x=0,y=Correlation)) +
  geom_segment(aes(x = 0, y = -0.05, xend = 0, yend = 1.05),
               arrow = arrow(length = unit(0.3, "cm"), end="both"), size=1, color="#1d47a1") +
  geom_segment(aes(x=0, y=Correlation, xend=1, yend=pheno_point), alpha=0.3) +
  geom_point(shape=1, size=3) +
  geom_point(aes(x=1,y=pheno_point),size=0.1) +
  scale_y_continuous(limits=c(-0.05,1.05), expand=c(0,0)) +
  scale_x_continuous(limits=c(-0.1,1.01)) +
  theme(plot.margin = margin(48,0,19,5), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_blank(), axis.title.y = element_text(size=12), axis.text.y = element_text(size=11)) +
  labs(x="")

lay <- rbind( c(1,2,2,2,2,2))
p <- grid.arrange(p2, p1, ncol = 2, layout_matrix=lay, 
                  top=textGrob("Genetic Correlation by Relative Heritability", gp=gpar(fontsize=16)))
p  


a <- c(1,2,3,4,5,6,7,8,9,10)
a[1:3+2]
