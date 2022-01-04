require("dplyr")
require("tidyr")
require("reshape2")
require("matrixStats")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

pheno <- "height"
setwd("~/Research/GWAS-frontera/OLD/GWAS_Results_OLD/height")
df_names <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t")

### simulation
setwd("~/Research/GWAS-frontera/mash/simulation")
get_df <- function(snps,h,e_ratio) {
  df <- read.csv(paste0(snps,"_",h,"_",e_ratio,".txt"), sep="\t")
  df <- data.frame(cbind(df_names$mix_0, df$x))
  colnames(df) <- c("Name", "Mean")
  df_values <- df
  
  # split matrice names
  df_values <- df_values %>%
    separate(Name, c("sex","correlation","effect"), sep="[_]", fill="right") %>%
    mutate(effect = paste0(sex, effect))
  return(df_values)
}

prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}

# split between null and values
snp_list <- c(100,1000,10000)
h <- "0.5"
e_ratio <- "1"

for (snp in snp_list) {
  df_values <- get_df(snp, h, e_ratio)
  ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4)])
  null <- prepare_df(df_values[1,c(2,3,4)])
  null$Mean <- as.numeric(null$Mean) ; ave$Mean <- as.numeric(ave$Mean)
  if (snp == snp_list[1]) {
    df_ave <- ave ; df_null <- null
  } else {
    df_ave <- cbind(df_ave, ave$Mean) 
    df_null <- rbind(df_null, null) 
  }
}

colnames(df_ave) <- c('correlation','effect',sprintf("snps_%s",snp_list))
df_ave <- df_ave %>% filter(snps_100+snps_1000+snps_10000 !=0) %>% 
  melt(id.vars=c('correlation','effect'))
df_null$parameter <- snp_list

big <- ggplot(df_ave, aes(x= effect, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(value,3)), color= "white", size= 4) +
  scale_y_continuous("correlation") +
  labs(title="Weights of Hypothesis Matrices", xlab="Magnitude", ylab="Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=12), legend.position = "none") +
  scale_fill_material("blue-grey") +
  facet_wrap(~variable,ncol=3)

small <- ggplot(df_null, aes(x= 0, y= 0, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(Mean,3)), color= "white", size= 4) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title="Weight of No Effect Matrice", xlab="Magnitude", ylab="Correlation") +
  theme_pubclean() +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position = "none",
        axis.ticks = element_blank(), plot.title = element_text(size=10)) +
  scale_fill_material("blue-grey") +
  facet_wrap(~parameter, ncol=3)
df_null  
small
lay <- rbind( c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(3,3,3,3))
p <- gridExtra::grid.arrange(big, small, ncol=1, layout_matrix=lay)
