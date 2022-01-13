require("dplyr")
require("tidyr")
require("reshape2")
require("matrixStats")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

# get covariance matrice names
pheno <- "height"
setwd("~/Research/GWAS-frontera/OLD/GWAS_Results_OLD/height")
df_names <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t")

### simulation
setwd("~/Research/GWAS-frontera/mash/simulation/ls6")

# get simulation mixture weight dataframe
get_df <- function(snps,h,e_ratio) {
  df <- read.csv(paste0(snps,"_",h,"_",e_ratio,".txt"), sep="\t")
  df <- data.frame(cbind(df_names$mix_0, df$x))   # combine names with dataframe
  colnames(df) <- c("Name", "Mean")
  df_values <- df
  
  # split matrice names
  df_values <- df_values %>%
    separate(Name, c("sex","correlation","effect"), sep="[_]", fill="right") %>%
    mutate(effect = paste0(sex, effect))
  return(df_values)
}

# order magnitude into factors, order by correlation and magnitude
prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}

##################################
# split between null and values
#snp_list <- c(100,1000,10000)
snp <- "100"
#h_list <- c("0.05","0.1","0.5")
h <- "0.5"                             
#e_ratio <- "1"   
e_ratio_list <- c("1","1.5","5")

#for (snp in snp_list) {
#for (h in h_list) {
for (e_ratio in e_ratio_list) {
  df_values <- get_df(snp, h, e_ratio)
  ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4)])
  null <- prepare_df(df_values[1,c(2,3,4)])
  null$Mean <- as.numeric(null$Mean) ; ave$Mean <- as.numeric(ave$Mean)
  #if (snp == snp_list[1]) {
  #if (h == h_list[1]) {
  if (e_ratio == e_ratio_list[1]) {
    df_ave <- ave ; df_null <- null
  } else {
    df_ave <- cbind(df_ave, ave$Mean) 
    df_null <- rbind(df_null, null) 
  }
}

# rename dataframe columns
#colnames(df_ave) <- c('correlation','effect',sprintf("snps_%s",snp_list))
#colnames(df_ave) <- c('correlation','effect',sprintf("h2_%s",h_list))
colnames(df_ave) <- c('correlation','effect',sprintf("E_ratio_%s",e_ratio_list))

# filter out weights that are 0 throughout snp numbers
df_ave <- df_ave %>% 
#  filter(snps_100+snps_1000+snps_10000 !=0) %>% 
  melt(id.vars=c('correlation','effect'))

#df_null$parameter <- snp_list
#df_null$parameter <- h_list
df_null$parameter <- e_ratio_list

big <- ggplot(df_ave, aes(x= effect, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(value,3)), color= "white", size= 4) +
  scale_y_continuous(breaks=seq(-1,1,0.25)) +
  theme_pubclean() +
  theme(axis.text=element_text(size=12), legend.position = "none", 
        plot.title = element_text(size=16), axis.title = element_text(size=14)) +
  labs(title="Weights of Hypothesis Matrices", x="Magnitude", y = "Correlation") +
  scale_fill_material("blue-grey") +
  facet_wrap(~variable,ncol=3)


small <- ggplot(df_null, aes(x= 0, y= 0, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(Mean,3)), color= "white", size= 4) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title="Weight of No Effect Matrice") +
  theme_pubclean() +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position = "none",
        axis.ticks = element_blank(), plot.title = element_text(size=10)) +
  scale_fill_material("blue-grey") +
  facet_wrap(~parameter, ncol=3) +
  labs(caption= paste0("# causal SNPs = ",snp,"\n",
                       "heritability = ", h, "\n"
                       #, "Environmental Variance Ratio = ",e_ratio
                       ))

lay <- rbind( c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(3,3,3,3))
p <- gridExtra::grid.arrange(big, small, ncol=1, layout_matrix=lay)


