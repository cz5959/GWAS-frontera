library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(ggsci)

setwd("~/Research/GWAS-frontera/PGS")

df <- read.csv("PGS_results_4.txt", sep="\t")
head(df)
# inc r2
df <- df[c(-3,-4,-7,-8,-11,-12,-15,-16)]
# r2
#df <- df[c(-5,-6,-9,-10,-13,-14)]
df$Sex[df$Sex == 'male'] <- 'm'
df$Sex[df$Sex == 'female'] <- 'f'

df_combined <- df[df$Sex == 'combined',]
df_sep <- df[df$Sex != 'combined',]
df_sep$Phenotype <- paste0(df_sep$Phenotype,"_",df_sep$Sex)

melt_out <- function(df) {
  df1 <- melt(df[c(1,2,3,5,7,9)], id.vars = c('Phenotype','Sex'),
              variable.name = "r2_var", value.name = "r2_val")
  df2 <- melt(df[c(1,2,4,6,8,10)], id.vars = c('Phenotype','Sex'),
              variable.name = "r2_var", value.name = "r2_se_val")
  df2$r2_var <- df1$r2_var
  df <- merge(df1,df2, by=c('Phenotype','Sex','r2_var'))
  df <- df[order(df$Phenotype, decreasing = FALSE),]
  return(df)
}
head(df_sep)
df_sep <- melt_out(df_sep)
df_combined <- melt_out(df_combined)
# order factor
df_combined$r2_var <- factor(df_combined$r2_var, levels=c("m_incr2","as_incr2","ah_incr2","ab_incr2"))
df_combined <- df_combined[order(df_combined$r2_var),]

# # formate df_combined for analysis
# df_combined = dcast(df_combined, Phenotype ~ r2_var, value.var="r2_val")
# #write.table(df_combined, file = "combined_r2.txt", sep="\t", row.names=FALSE, quote = FALSE)
# # ttest for difference between additive and mash
# df_add <- df_combined[df_combined$r2_var == "_r2", c(1,2)]
# df_mash <- df_combined[df_combined$r2_var == "m_r2",c(1,4)]
# df_add <- df_add %>%
#   mutate(mean_diff = r2_val - df_mash$r2_val) %>%
#   mutate(se_diff = sqrt((df_add$r2_se_val)^2 + (df_mash$r2_se_val)^2) ) %>%
#   mutate(ttest = mean_diff - se_diff)
# #write.table(df_add, file = "combined_diff.text", sep="\t", row.names=FALSE, quote = FALSE)

# combined plot, r2
#pdf(file="pgs_comparison_five_combined.pdf",width=8,height=14)
p <- ggplot(data=df_combined, aes(x=Phenotype, y=r2_val, fill=r2_var)) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-r2_se_val, ymax=r2_val+r2_se_val), alpha= 0.8, show.legend = FALSE,
                position='dodge', stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  labs(title='PGS Comparison Over Five Folds - Combined Models', y="R2", fill="Model") +
  coord_flip() + 
  theme_classic() + 
  scale_fill_npg(labels = c('additive both-sex', 'additive same-sex', 'mash',"additive both-sex\n half sample size")) +
  theme(axis.text = element_text(size=9), axis.title = element_text(size=11), plot.title=element_text(size=14),
        legend.title=element_text(size=11), legend.text=element_text(size=9), 
        legend.position = 'top', legend.background = element_rect(linetype='solid', color='black'))
p
#dev.off()

## combined plot, edited 3SE
ggplot(data=df_combined, aes(x=Phenotype, y=r2_val, fill=(r2_var))) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-3*r2_se_val, ymax=r2_val+3*r2_se_val), alpha= 0.8, 
                show.legend = FALSE, position='dodge', stat='identity', size=0.2) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title='PGS Comparison', y="Incremental R2", fill="Model") +
  coord_flip() + 
  theme_classic() + 
  scale_fill_manual(values = c("#2b62d9", "#d1b724", "#649e98","#b0464f")) +
  theme(axis.text = element_text(size=9), axis.title.x = element_text(size=11),
        plot.title=element_text(size=14), axis.title.y = element_blank(),
        legend.position = "none") +
  annotate("text", x = 26, y=0.105, label = "Additive", hjust = 1, color="#b0464f", size=3.4 ) +
  annotate("text", x = 25, y=0.105, label = "Additive (half sample)", hjust = 1, color="#649e98", size=3.4 ) +
  annotate("text", x = 24, y=0.105, label = "Sex-specific additive", hjust = 1, color="#d1b724", size=3.4) +
  annotate("text", x = 23, y=0.105, label = "Sex-specific covariance aware", hjust = 1, color="#2b62d9", size=3.4 )
  
  



############### SAMPLE SIZE

setwd("~/Research/GWAS-frontera/PGS")
df <- read.csv("PGS_samplesize.txt", sep="\t")
df2 <- read.csv("PGS_samplesize_intercept.txt", sep="\t")
df2$perc <- 50
head(df2)

ggplot(df, aes(x=perc, y=ah_inc_r2, color=phenotype)) +
  geom_point() +
  geom_errorbar(aes(ymin=ah_inc_r2-3*ah_inc_r2_se, ymax=ah_inc_r2+3*ah_inc_r2_se),
                width=0.2, position=position_dodge(0.05)) +
  geom_rect(aes(xmin=48, xmax=102, ymin=ah_inc_r2-3*ah_inc_r2_se, 
                ymax=ah_inc_r2+3*ah_inc_r2_se, fill=phenotype), 
            data=df2, alpha=0.2, linetype=0, show.legend = FALSE) +
  geom_hline(aes(yintercept=ah_inc_r2, color=phenotype, linetype=model), df2) +
  labs(x="% total sample size", y="incremental R^2") +
  theme_classic() +
  scale_color_manual(values = c("#b0464f","#2b62d9")) +
  scale_fill_manual(values = c("#b0464f","#2b62d9")) +
  scale_linetype_manual(values=c("longdash","dotted"), 
                        labels=c("Sex-specific additive","Sex-specific covariance aware")) +
  scale_x_continuous(expand=c(0,0), limits=c(48,102), breaks=seq(50,100,10)) +
  theme(legend.position="top") +
  guides(linetype=guide_legend(nrow=2), color=guide_legend(nrow=2))

  

