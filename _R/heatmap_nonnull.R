require("dplyr")
require("tidyr")
require("reshape2")
require("matrixStats")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

setwd("~/Documents/Harpak/GxSex/mash/")
#setwd("~/Research/GWAS-frontera/mash/")
m_names <- read.csv("matrice_names.txt", sep="\t")

setwd("~/Documents/Harpak/GxSex/mash/simulation/nonnull2")
#setwd("~/Research/GWAS-frontera/mash/simulation/nonnull2")
# f=18
# m=64
# s=2
# snps=1000
# same="same"
# 
# name <- paste0("mash_",snps,"_nonnull_",f,"_",m,"_",s,"_",same)
name <- "mash_1000_0_all"
df <- read.csv(paste0(name,".txt"), sep="\t")
mean <- rowMeans(df)
se <- rowSds(as.matrix(df)) / sqrt(length(colnames(df)) - 1)

df <- data.frame(cbind(m_names[1:4], mean, se))
df$effect <- paste0(df$sex, df$magnitude)
prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}
df_nonnull <- prepare_df(df[2:nrow(df),c(3,5,6,7)])
df_null <- prepare_df(df[1,c(3,4,5,6,7)])

# format labels
df_nonnull <- df_nonnull %>% 
  mutate(mean = as.numeric(mean)) %>%
  mutate(mean_lab = ifelse(mean < 0.0005, "0%", sprintf("%.1f%%", round(mean*100,1)) )) %>%
  mutate(se = as.numeric(se)) %>%
  mutate(se_lab = ifelse(se < 0.0005, "0%", sprintf("%.1f%%", round(se*100,1)) )) 
df_null <- df_null %>% 
  mutate(mean = as.numeric(mean)) %>%
  mutate(mean_lab = ifelse(mean < 0.0005, "0%", sprintf("%.1f%%", round(mean*100,1)) )) %>%
  mutate(se = as.numeric(se)) %>%
  mutate(se_lab = ifelse(se < 0.0005, "0%", sprintf("%.1f%%", round(se*100,1)) )) 


### BIG PLOT
effect_labels <-  c('female-\nspecific','female x3', 'female x2', 'female x1.5','equal','male x1.5','male x2','male x3','male-\nspecific')

#png(file=paste0(name,".png"), width=6.5, height=4.8, units="in", res=300)
pdf(file=paste0(name,".pdf"), width=6.5, height=4.8)
big <- ggplot(df_nonnull, aes(x= effect, y= correlation, fill= mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=mean_lab), color= "black", size= 2.7, vjust=-0.1) +
  geom_text(aes(label=paste("\u00B1",se_lab)), color= "black", size= 2.2, vjust=1.5) +
  scale_y_continuous(breaks=seq(-1,1,0.25), expand=c(0,0)) +
  scale_x_discrete(labels= effect_labels) +
  #labs(title=paste0(f,"% female; ",m,"% male; x",s,"; snps=",snps,"; same sign")) +
  xlab("Magnitude") + ylab("Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=9), axis.title = element_text(size=11), plot.title = element_text(size=14), 
       legend.position = "none") +
  scale_fill_gradient(low="gray98",high="#829ed9")

small <- ggplot(df_null, aes(x= 0, y= 0, fill= mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=mean_lab), color= "white", size= 2.7, vjust=-0.1) +
  geom_text(aes(label=paste("\u00B1",se_lab)), color= "white", size= 2.2, vjust=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Weight of \nNo Effect Matrix") +
  theme_pubclean() +
  theme(axis.text=element_blank(), axis.title=element_blank(), 
        axis.title.y = element_text(size=10, angle=360, vjust=0.5),
        legend.position = "none",axis.ticks = element_blank(), plot.title = element_blank()) +
  scale_fill_gradient(low="#2b62d9",high="#2b62d9")

lay <- rbind( c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(2,2,3,3,3))
p <- gridExtra::grid.arrange(big, small, ncol=1, layout_matrix=lay)


dev.off()

### SMALL PLOT
nan_weight <- 1 / (1 - df_null$mean[1])
df_nonnull$mean = df_nonnull$mean * nan_weight
df_nonnull$effect <- as.character(df_nonnull$effect) ; df_nonnull$correlation <- as.numeric(df_nonnull$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_nonnull %>% filter(substr(effect,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(mean)) %>%
    as.data.frame()
  return(df_sex)
}
for (s in c('f','m','e')) {
  assign(s, group_sex(s))
}

df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2])); colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
  c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
  c('partial,\npositive', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
  c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
  c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))
colnames(df_small) <- c('corr', 'female > male', 'female = male', 'female < male')
df_small <- melt(df_small, id.vars=c('corr'))

df_small$corr <- factor(df_small$corr, levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))
df_small$variable <- factor(df_small$variable, levels = c('female > male', 'female = male', 'female < male'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)

colnames(df_small) <- c("correlation", "magnitude", "value")

sum_mag <- df_small %>%
  group_by(magnitude) %>%
  summarise(sum= sum(value))
sum_corr <- df_small %>%
  group_by(correlation) %>%
  summarise(sum = sum(value))

# PLOT 
#png(file=paste0(name,"_small.png"), width=3.5, height=2.3, units="in", res=200)
pdf(file=paste0(name,"_small.pdf"), width=3.5, height=2.3)
ggplot(df_small, aes(x=magnitude, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=sprintf("%.0f%%", (value*100))), color= "black", size= 2.5) +
  labs(title="Covariance of Genetic Effects: Compact Representation") +
  xlab("Magnitude") + ylab("Correlation") +
  theme_pubclean() +
  # theme(axis.text=element_text(size=9), axis.title = element_text(size=10, hjust=0.5),
  #        plot.title = element_text(size=11, hjust=1), legend.position = "none") +
  theme(axis.text=element_text(size=8), axis.title = element_text(size=9, hjust=0.5),
        plot.title = element_blank(), legend.position = "none") +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low="gray98",high="#829ed9") + 
  annotate("text", x=1:3, y = 0.4, size = 2.8, label = sprintf("%.0f%%", (sum_mag$sum*100))) +
  annotate("text", x = 3.65, y=1:4, size = 2.8, label = sprintf("%.0f%%", (sum_corr$sum*100))) +
  coord_cartesian(xlim=c(1,3.15), ylim=c(0.9,3.9), clip="off")

dev.off()


