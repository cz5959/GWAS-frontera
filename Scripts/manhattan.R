#!/usr/bin/env Rscript

libLoc <- "/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/"
library("optparse", lib.loc=libLoc)
library("crayon",lib.loc=libLoc)
library("dplyr",lib.loc=libLoc)
library("tidyr",lib.loc=libLoc)
library("withr",lib.loc=libLoc)
library("ggplot2",lib.loc=libLoc)
library("grid",lib.loc=libLoc)
library("gridExtra",lib.loc=libLoc)
library("labeling",lib.loc=libLoc)
library("farver",lib.loc=libLoc)
library("digest",lib.loc=libLoc)
library("backports",lib.loc=libLoc)
library("ggpubr",lib.loc=libLoc)
library("ggsci",lib.loc=libLoc)

# user input phenotype
option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
  make_option(c("-n","--name"), type="character", default=NULL,help="formatted phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)
pheno_name <- opt$name

setwd( paste0("/scratch1/08005/cz5959/GWAS_Results/",pheno) )
male_name <- paste0("male_all.",pheno,".glm.linear")
female_name <- paste0("female_all.",pheno,".glm.linear")
male_df <- read.csv(male_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ))
female_df <- read.csv(female_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ))

edit_df <- function(df) {
  for (p in c(5e-30, 5e-20, 5e-10, 5e-5, 0.05, 0.1, 0.5)) {
    df <- bind_rows(filter(df, P<=p), sample_frac(filter(df, P>p), 0.3))
  }
  plot_df <- df %>% 
  drop_na() %>%
  group_by(X.CHROM) %>% 
  summarize(CHR_LEN=max(POS)) %>%
  mutate(TOT=cumsum(as.numeric(CHR_LEN))-CHR_LEN) %>%
  select(-CHR_LEN) %>%
  left_join(df, ., by=c("X.CHROM"="X.CHROM")) %>%
  arrange(X.CHROM, POS) %>%
  mutate(POS_CUM=POS+TOT) %>%
  mutate(COLOR=ifelse(X.CHROM %% 2, 1,2))
  return(plot_df)
}

x_axis <- function(plot_df) {
  axis_df <- plot_df %>% 
    group_by(X.CHROM) %>%
    summarize(CENTER=(max(POS_CUM)+min(POS_CUM))/2 )
  return(axis_df)
}

female_df <- edit_df(female_df)
male_df <- edit_df(male_df)
axis_df <- x_axis(female_df)

pdf(file=paste0(pheno,"_miami.pdf"), width=10, height=4)

female_plot <- ggplot(female_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=0.1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_continuous(expand=c(0,0)) +
  
  theme_pubclean() +
  theme(legend.position = "none", axis.title = element_blank(), plot.margin = margin(20,20,0,20),
        axis.text = element_text(size=10), panel.grid.major.y=element_blank(), plot.title=element_text(size=16, hjust=0.5)) +
  scale_color_manual(values=c("#d67629","#1d47a1")) +
  labs(title=pheno_name)

male_plot <- ggplot(male_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=0.1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_reverse(expand=c(0,0)) +
  
  theme_pubclean() + 
  theme(legend.position = "none",  axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
        plot.margin =  margin(0,20,20,20), panel.grid.major.y=element_blank(), 
        axis.text.y = element_text(size=10), axis.title.y = element_blank(), axis.title.x = element_text(size=12)) +
  scale_color_manual(values=c("#d67629","#1d47a1")) +
  labs(x="Chromosome")

p <- grid.arrange(female_plot, male_plot, nrow = 2, left=textGrob(bquote(-log[10] ~P), rot=90, vjust=1, gp=gpar(fontsize=12)))
p

dev.off()

