require("dplyr")
require("tidyr")
require("ggplot2")
require("grid")
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
  theme(legend.position = "none", axis.title = element_blank(), plot.margin = margin(20,5.5,0,5.5), panel.grid.major.y=element_blank(),
        axis.text = element_text(size=9), plot.title=element_text(size=16, hjust=0.5)) +
  scale_color_manual(values=c("#d67629","#d6a88c")) +
  labs(title=pheno_name) +
  annotation_custom(grobTree(textGrob("female", x=0.8, y=0.8, gp = gpar(col="#d67629", fontsize=11) )))

male_plot <- ggplot(male_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=1) +
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_reverse(expand=c(0,0)) +
  theme_pubclean() + 
  theme(legend.position = "none",  axis.text.x=element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = margin(0,5.5,40,5.5), panel.grid.major.y=element_blank(),
        axis.text.y = element_text(size=9), axis.title.x = element_text(size=11)) +
  scale_color_manual(values=c("#207335","#99c4a3")) +
  labs(x="Chromosome") +
  annotation_custom(grobTree(textGrob("male", x=0.8, y=0.2, gp = gpar(col="#207335", fontsize=11) )))

p <- grid.arrange(female_plot, male_plot, nrow = 2, left=textGrob( bquote(-log[10] ~P) ,
                                                                  rot=90, vjust=1, gp=gpar(fontsize=11)))

#dev.off()



