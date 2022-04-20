require("dplyr")
require("tidyr")
require("tidyverse")
require("reshape2")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

# load mixture names
setwd("~/Research/GWAS-frontera/mash/mash_mixprop_100")
mix <- read.csv("heightmixprop_100_all.txt", sep="\t")
mix <- mix$mix_0

# phenotype list
setwd("~/Research/GWAS-frontera/mash pvalue")
pheno <- "testosterone"; title <- "Testosterone"
p_values <- factor(c("5e-08", "1e-05", "0.05", "1"))

# initialize list for null weight
#all_null <- data.frame(pval=p_values)

null <- NULL  # START LOOP

for ( p_val in p_values ) {

# load mixture proportions by p-value
df <- read.csv(paste0(pheno,"_",p_val,"_same.txt"), sep="\t")
df <- cbind(mix, df)

# ORGANIZE
# split matrice names
df <- df %>%
  separate(mix, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
  mutate(magnitude = paste0(sex, magnitude))

# organize by correlation and magnitude
prepare_df <- function(df) {
  df$magnitude <- factor(df$magnitude, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, magnitude)
  return(df)
}

# split between null and values
df_ave <- prepare_df(df[2:nrow(df),c(2,3,4)])
df_null <- prepare_df(df[1,c(2,3,4)])
null <- append(null, df_null$x)

# CONDENSE
nan_weight <- 1 / (1 - df_null$x[1])
df_ave$x = df_ave$x * nan_weight
df_ave$magnitude <- as.character(df_ave$magnitude) ; df_ave$correlation <- as.numeric(df_ave$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_ave %>% filter(substr(magnitude,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(x)) %>%
    as.data.frame()
  return(df_sex)
}
for (s in c('f','m','e')) {
  assign(s, group_sex(s))
}

# group by correlation
df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2])); colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
  c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
  c('partial,\npositive', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
  c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
  c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))

# edit col names and data type
colnames(df_small) <- c('corr', 'female > male', 'female = male', 'female < male')
df_small <- melt(df_small, id.vars=c('corr'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)
colnames(df_small) <- c("correlation", "magnitude", "value")

# sum across correlation and magnitude groups
sum_mag <- df_small %>%
  group_by(magnitude) %>%
  summarise(sum= sum(value))
sum_corr <- df_small %>%
  group_by(correlation) %>%
  summarise(sum = sum(value))

# assign df to variable
assign(paste0(p_val,"_mag"), as.data.frame(sum_mag))
assign(paste0(p_val,"_corr"), as.data.frame(sum_corr))

}
## END OF LOOP

# no effect weight

all_null <- cbind(all_null, null)
all_null
colnames(all_null) <- c("pval", "Height","BMI","Waist:hip (bmi adjusted)","Testosterone")
all_null <- all_null %>%
  gather(key=Phenotype, value=Weight, 2:5) %>%
  mutate(pval = factor(pval, levels = c("5e-08", "1e-05", "0.05", "1")))
#write.table(all_null, file="noeffect_weight_same.txt",sep="\t",row.names = F, quote=F)

# NO EFFECT PLOT
# 3 x 6
colors <- c("#d1b724", "#563f61", "#b0464f", "#2b62d9")
ggplot(all_null, aes(x=pval, y=Weight, group=Phenotype, color=Phenotype)) +
  geom_point() +
  geom_line() +
  labs(x="P-value Threshold") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(axis.title=element_text(size=10), axis.text=element_text(size=9))

# ORGANIZE FOR PLOT
# merge p-value mixture weight tables
corr_list <- list(`0.05_corr`, `1_corr`,`1e-05_corr`,`5e-08_corr`)
corr_all <- corr_list %>% 
  reduce(full_join, by='correlation') %>%
  rename("5e-8"=sum.y.y, "1e-5"=sum.x.x, "0.05"=sum.x, "1"=sum.y) %>%
  gather(threshold, weight, 2:5) %>%
  mutate(threshold = factor(threshold, levels=c("5e-8", "1e-5", "0.05", "1")))
  
corr_all$correlation <- factor(corr_all$correlation, 
                               levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))

mag_list <- list(`0.05_mag`, `1_mag`,`1e-05_mag`,`5e-08_mag`)
mag_all <- mag_list %>% 
  reduce(full_join, by='magnitude') %>%
  rename("5e-8"=sum.y.y, "1e-5"=sum.x.x, "0.05"=sum.x, "1"=sum.y) %>%
  gather(threshold, weight, 2:5) %>%
  mutate(threshold = factor(threshold, levels=c("5e-8", "1e-5", "0.05", "1")))


# PLOT
# pdf 3x4
colors <- c("#d1b724", "#563f61", "#b0464f", "#2b62d9")
#pdf(file=paste0(pheno,"_corr_same.pdf"), width=4, height=3)
corr_plot <- 
  ggplot(corr_all, aes(x=threshold, y=weight*100, fill=correlation)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(title=title, x="P-value Threshold", y="Weight", fill="Correlation") +
  scale_fill_manual(values = colors) +
  theme(axis.title=element_text(size=11), axis.text=element_text(size=9))
print(corr_plot)
#dev.off()

colors <- c("#d67629", "#829ed9", "#207335")
pdf(file=paste0(pheno,"_mag_same.pdf"), width=4, height=3)
mag_plot <- 
  ggplot(mag_all, aes(x=threshold, y=weight*100, fill=magnitude)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(title=title, x="P-value Threshold", y="Weight", fill="Magnitude") +
  scale_fill_manual(values = colors) +
  theme(axis.title=element_text(size=11), axis.text=element_text(size=9))
print(mag_plot)
dev.off()





