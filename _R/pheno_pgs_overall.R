

setwd("~/Research/Phenotypes/")
pheno_names <- read.csv("pheno_names.txt", sep="\t")

setwd("~/Research/GWAS-frontera/PGS")
df <- read.csv("sexspecific_pheno_pgs_lm.txt", sep="\t")
head(df)
df <- merge(pheno_names, df, by.x="Code", by.y="Phenotype")
df_m <- df[df$PGS_sex == "M",]
df_f <- df[df$PGS_sex == "F",]

#
odd <- seq(1,54,2); even <- seq(2,54,2)
df_m <- df_m %>%
  mutate(Beta_ratio = rep(abs(Beta[even] / Beta[odd]), each=2)) %>%
  mutate(Beta_ratio_err = rep((Beta_ratio[odd])*sqrt((Err[odd]^2/Beta[odd]^2)+
                                                       (Err[even]^2/Beta[even]^2)), each=2) ) %>%
  mutate(outlier = ifelse(Phenotype == "Testosterone", 1, 0))

df_f <- df_f %>%
  mutate(Beta_ratio = rep(abs(Beta[even] / Beta[odd]), each=2)) %>%
  mutate(Beta_ratio_err = rep((Beta_ratio[odd])*sqrt((Err[odd]^2/Beta[odd]^2)+
                                                       (Err[even]^2/Beta[even]^2)), each=2) ) %>%
  mutate(outlier = ifelse(Phenotype == "Testosterone", 1, 0))


df_f <- slice(df_f, odd)
df_m <- slice(df_m, odd)
pheno_order <- reorder(df_f$Phenotype, df_f$Beta_ratio)

### PLOT ###
# MALE PGS 4x5.5
df_m1 <- df_m[df_m$Code != "testosterone",] 
df_m2 <- df_m[df_m$Code == "testosterone",] 
rects <- data.frame(xstart = seq(0.5,27.5,1), xend = seq(1.5,28.5,1))
rects <- rects[1:27,]
g1 <- ggplot(df_m1, aes(x=pheno_order[c(1:19,21:27)], y=Beta_ratio)) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=Beta_ratio-(2*Beta_ratio_err), ymax=Beta_ratio+(2*Beta_ratio_err)), 
                alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(size=9),
        legend.position = "none") +
  scale_y_continuous(expand=c(0,0)) + 
  geom_rect(data=rects, aes(xmin=xstart-1, xmax=xend-1, ymin=0.6, ymax=2.1), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  coord_flip()

g2 <- ggplot(df_m2, aes(x=Phenotype, y=Beta_ratio), color="#2b62d9") +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=Beta_ratio-(2*Beta_ratio_err), ymax=Beta_ratio+(2*Beta_ratio_err)), 
                alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text.y= element_blank(), axis.text.x=element_text(size=9),
        legend.position = "none", axis.line.y=element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(expand=c(0,0), breaks=c(40,100,160)) +
  geom_rect(data=rects, aes(xmin=xstart, xmax=xend, ymin=10, ymax=200), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  coord_flip()

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g2)
maxHeight = grid::unit.pmax(gA$heights[2:5], gB$heights[2:5])
gA$heights[2:5] <- as.list(maxHeight) ; gB$heights[2:5] <- as.list(maxHeight)
plot <- grid.arrange(gA, gB, ncol=2, nrow=1, widths=c(5,1))



# FEMALE PGS 4x5.5
rects <- data.frame(xstart = seq(0.5,27.5,1), xend = seq(1.5,28.5,1))
rects <- rects[1:27,]
ggplot(df_f, aes(x=pheno_order, y=Beta_ratio)) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=Beta_ratio-(2*Beta_ratio_err), ymax=Beta_ratio+(2*Beta_ratio_err)), 
                alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(size=9),
        legend.position = "none") +
  scale_y_continuous(expand=c(0,0), breaks = seq(-0.5,2,0.5)) +
  #scale_y_reverse(expand=c(0,0)) +
  #scale_x_discrete(position="top") +
  geom_rect(data=rects, aes(xmin=xstart, xmax=xend, ymin=-0.9, ymax=2), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  coord_flip()

 
### STAT TEST ###
## spearman correlation between panels
cor.test(df_m$Beta_ratio, df_f$Beta_ratio, method="spearman")

#
df_m <- df_m %>%
  mutate(is_sig = ifelse((Beta_ratio + 1.645*Beta_ratio_err < 1) | (Beta_ratio - 1.645*Beta_ratio_err > 1),1,0 ))
df_f <- df_f %>%
  mutate(is_sig = ifelse((Beta_ratio + 1.645*Beta_ratio_err < 1) | (Beta_ratio - 1.645*Beta_ratio_err > 1),1,0 ))

# number traits concordantly significant: 11
# traits both significant but not concordant: 1; waist circumference
# both waist to hip ratios opposite signs in mean Beta ratio

## two sample t-test
df_mm <- df_m[df_m$Sex == 1,]; df_mf <- df_m[df_m$Sex == 0,]
df_fm <- df_f[df_f$Sex == 1,]; df_ff <- df_f[df_f$Sex == 0,]
# male polygenic score
df_mf$t <- (df_mm$Beta - df_mf$Beta) / sqrt( df_mm$Err^2 + df_mf$Err^2)
df_mf$p <- pt(abs(df_mf$t), 99998, lower.tail = F)
(df_mf[df_mf$p < 0.05, c("Code")])
# female polygenic score
df_ff$t <- (df_fm$Beta - df_ff$Beta) / sqrt( df_fm$Err^2 + df_ff$Err^2)
df_ff$p <- 2*pt(abs(df_ff$t), 99998, lower.tail = F)
(df_ff[df_ff$p < 0.05, c("Code")])

