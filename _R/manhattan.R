require("dplyr")
require("tidyr")
require("ggplot2")
require("gridExtra")
require("ggpubr")

# https://www.r-graph-gallery.com/101_Manhattan_plot.html
# https://github.com/juliedwhite/miamiplot/blob/master/vignettes/scratch_miamiplots.Rmd

pheno <- 'height'
pheno_name <- "Height"
setwd("~/Research/GWAS-frontera/GWAS_Results/height")
file_name <- paste0("male_sample.",pheno,".glm.linear")

male_df <- read.csv(file_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ),nrow=1000)
female_df <- read.csv(file_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ),nrow=1000)

edit_df <- function(df) {
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

#pdf(file=paste0(pheno,"_miami.pdf"), width=8, height=4)

female_plot <- ggplot(female_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_continuous(expand=c(0,0)) +
  
  theme_pubclean() +
  theme(legend.position = "none", axis.title.x = element_blank(), plot.margin = margin(20,20,0,20),
        axis.text = element_text(size=11), panel.grid.major.y=element_blank()) +
  scale_color_npg() +
  labs(title=paste0("Miami Plot - ", pheno_name), y="Female")

male_plot <- ggplot(male_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_reverse(expand=c(0,0)) +
  
  theme_pubclean() + 
  theme(legend.position = "none",  axis.text.x=element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(0,20,20,20),
        axis.text = element_text(size=11), panel.grid.major.y=element_blank()) +
  scale_color_npg() +
  labs(x="Chromosome",y="Male")

p <- grid.arrange(female_plot, male_plot, nrow = 2, left=textGrob("-log10(P)",rot=90,vjust=1))
p
#dev.off()




