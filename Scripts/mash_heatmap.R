require("dplyr")
require("tidyr")
require("reshape2")
require("matrixStats")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

pheno <- "testosterone"
setwd(paste0("~/Research/GWAS-frontera/OLD/GWAS_Results_OLD/",pheno))

df <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t")

# get mean and sem (NOT SIM)
df_values <- data.frame( Name = df$mix_0, Mean = rowMeans(df[2:101]), 
                         SE = rowSds(as.matrix(df[2:101])) / sqrt(length(colnames(df)) - 1) )
# split matrice names
df_values <- df_values %>%
  separate(Name, c("sex","correlation","effect"), sep="[_]", fill="right") %>%
  mutate(effect = paste0(sex, effect))

prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}

# split between null and values
df_ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4,5)])
df_null <- prepare_df(df_values[1,c(2,3,4,5)])
df_null$Mean <- as.numeric(df_null$Mean)
df_ave$Mean <- as.numeric(df_ave$Mean)

effect_labels <-  c('female','female x 3', 'female x 2', 'female x 1.5','equal','male x 1.5','male x 2','male x 3','male')
big <- ggplot(df_ave, aes(x= effect, y= correlation, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(Mean,3)), color= "white", size= 6, vjust=-0.1) +
  geom_text(aes(label=round(SE,3)), color= "white", size= 4, vjust=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(labels= effect_labels) +
  labs(title="Weights of Hypothesis Matrices", xlab="Magnitude", ylab="Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=14), legend.position = "none") +
  scale_fill_material("blue-grey")

small <- ggplot(df_null, aes(x= 0, y= 0, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(Mean,3)), color= "white", size= 6, vjust=-0.1) +
  geom_text(aes(label=round(SE,3)), color= "white", size= 4, vjust=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title="Weight of No Effect Matrice", xlab="Magnitude", ylab="Correlation") +
  theme_pubclean() +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position = "none",
        axis.ticks = element_blank(), plot.title = element_text(size=12)) +
  scale_fill_material("blue-grey")

lay <- rbind( c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(2,3,3,3))
p <- gridExtra::grid.arrange(big, small, ncol=1, layout_matrix=lay)


##########################################################################
### SMALL HEATMAP ###
nan_weight <- 1 / (1 - df_null$Mean[1])
df_ave$Mean = df_ave$Mean * nan_weight
df_ave$effect <- as.character(df_ave$effect) ; df_ave$correlation <- as.numeric(df_ave$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_ave %>% filter(substr(effect,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(Mean)) %>%
    as.data.frame()
  return(df_sex)
}
for (s in c('f','m','e')) {
  assign(s, group_sex(s))
}

df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2])); colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
  c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
  c('partial', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
  c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
  c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))
colnames(df_small) <- c('corr', 'F > M', 'F = M', 'M > F')
df_small <- melt(df_small, id.vars=c('corr'))

df_small$corr <- factor(df_small$corr, levels = c('negative', 'uncorrelated', 'partial', 'perfect'))
df_small$variable <- factor(df_small$variable, levels = c('F > M', 'F = M', 'M > F'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)

ggplot(df_small, aes(x=variable, y= corr, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(value,3)), color= "white", size= 4) +
  labs(title="Weights of Hypothesis Matrices", xlab="Magnitude", ylab="Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=12), legend.position = "none") +
  scale_fill_material("blue-grey")



