library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(ggsci)

setwd("~/Research/GWAS-frontera/GWAS_Results")

df <- read.csv("pgs_linear_results_five.txt", sep="\t")
df <- df[c(-5,-6,-9,-10,-13,-14)]
df$Sex[df$Sex == 'male'] <- 'm'
df$Sex[df$Sex == 'female'] <- 'f'

df_combined <- df[df$Sex == 'combined',]
df_sep <- df[df$Sex != 'combined',]
df_sep$Phenotype <- paste0(df_sep$Phenotype,"_",df_sep$Sex)

melt_out <- function(df) {
  df1 <- melt(df[c(1,2,3,5,7)], id.vars = c('Phenotype','Sex'),
             variable.name = "r2_var", value.name = "r2_val")
  df2 <- melt(df[c(1,2,4,6,8)], id.vars = c('Phenotype','Sex'),
             variable.name = "r2_var", value.name = "r2_se_val")
  df2$r2_var <- df1$r2_var
  df <- merge(df1,df2, by=c('Phenotype','Sex','r2_var'))
  return(df)
}

df_sep <- melt_out(df_sep)
df_combined <- melt_out(df_combined)

# separated plot, r2
p <- ggplot(data=df_sep, aes(x=Phenotype, y=r2_val, fill=factor(r2_var))) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-r2_se_val, ymax=r2_val+r2_se_val), alpha= 0.8, show.legend = FALSE,
                position='dodge', stat='identity') +
  labs(title='PGS Comparison Over Five Folds', y="R2", fill="Model") + 
  coord_flip() + theme_classic() + scale_fill_npg(labels = c('additive both-sex', 'additive same-sex', 'mash')) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=14), plot.title=element_text(size=16),
        legend.title=element_text(size=12), legend.text=element_text(size=10), 
        legend.position = 'top', legend.background = element_rect(linetype='solid', color='black'))

p

# combined plot, r2
p <- ggplot(data=df_combined, aes(x=Phenotype, y=r2_val, fill=factor(r2_var))) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-r2_se_val, ymax=r2_val+r2_se_val), alpha= 0.8, show.legend = FALSE,
                position='dodge', stat='identity') +
  labs(title='PGS Comparison Over Five Folds - Combined Models', y="R2", fill="Model") + 
  coord_flip() + theme_classic() + scale_fill_npg(labels = c('additive both-sex', 'additive same-sex', 'mash')) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=14), plot.title=element_text(size=16),
        legend.title=element_text(size=12), legend.text=element_text(size=10), 
        legend.position = 'top', legend.background = element_rect(linetype='solid', color='black'))
p

