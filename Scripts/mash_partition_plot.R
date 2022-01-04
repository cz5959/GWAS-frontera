require("dplyr")
require("tidyr")
require("ggplot2")

setwd("~/Research/GWAS-frontera/mash/celltype")

### overlaps and cell type plot
df <- read.csv("celltype_sums.txt")
hist(df$x, breaks=seq(0.5,10.5,1), xaxp=c(1,11,10), xlim=c(0.5,10.5),
     xlab='# of Overlaps', main='Histogram of # of Overlaps')
df <- read.csv("celltype_colsums.txt", sep="\t")
barplot(df$x, names.arg=c(1:10), xlab='Cell Type', ylab='# SNPs',
        main='Number of SNPs per Cell Type')

### mash posterior weights distribution
df <- read.csv("testosterone_mash_weights_col24.txt", sep="\t")
head(df)
hist(df$m_0_1, 
     xlab='m_-m_0_1', main='Histogram of Testosterone Posterior Weights')

########## 
pheno <- "height"
df <- read.csv(paste0(pheno,"_mash_partitioned.txt"), sep="\t")

# split matrice names
df <- as.data.frame(t(df))
df <- df %>%
  mutate(matrice = rownames(df)) %>%
  separate(matrice, c("sex","correlation","effect"), sep="[_]", fill="right") %>%
  select(-c(13)) %>% 
  mutate(correlation = ifelse(substr(correlation,1,1) == '.',
                            paste0('-',substr(correlation,2,5)),correlation))
# group by sex
group_sex <- function(sex){
  df_sex <- df %>% 
    filter(substr(sex,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(across(V1:V10, mean)) %>%
    as.data.frame()
  return(df_sex)
}
# group by correlation
group_corr <- function(df_sex) {
  corr_df <- data.frame(rbind(
    c('perfect', colSums(df_sex[df_sex$correlation == 1, 2:11])),
    c('partial', colSums(df_sex[(df_sex$correlation > 0) & (df_sex$correlation < 1), 2:11])),
    c('uncorrelated', colSums(df_sex[(df_sex$correlation == 0), 2:11])),
    c('negative', colSums(df_sex[(df_sex$correlation < 0), 2:11]))
  ))
  return(corr_df)
}

for (s in c('f','m','e')) {
  assign(s, group_sex(s))
  assign(s, group_corr(get(s)))
}

celltypes <- c(" Adrenal Pancreas", "Cardiovascular", "CNS", "Connective Bone", "GI","Hematopoietic",
               "Kidney", "Liver","Other","Skeletal Muscle")
df_plot <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("corr","variable", "value", "celltype"))))
for (i in 2:11) {
  # make list of 10 cell type dataframes
  df_grouped <- data.frame(cbind(f[,1], f[,i], e[,i], m[,i]))
  colnames(df_grouped) <- c('corr', 'F > M', 'F = M', 'M > F')
  df_grouped <- melt(df_grouped, id.vars=c('corr'))
  df_grouped$corr <- factor(df_grouped$corr, levels = c('negative', 'uncorrelated', 'partial', 'perfect'))
  df_grouped$variable <- factor(df_grouped$variable, levels = c('F > M', 'F = M', 'M > F'))
  df_grouped <- df_grouped %>% 
    mutate_at(3, as.numeric) %>%
    arrange(corr, variable) %>%
    mutate(celltype=celltypes[i-1])
  df_plot <- rbind(df_plot, df_grouped)
}

# make graph
ggplot(df_plot, aes(x=variable, y= corr, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=round(value,3)), color= "white", size= 4) +
  theme_pubclean() +
  labs(title="Weights of Hypothesis Matrices") + ylab("Correlation") + xlab("Magnitude") +
  theme(axis.text=element_text(size=12), legend.position = "none") +
  scale_fill_material("blue-grey") +
  facet_wrap(~celltype,ncol=5)



