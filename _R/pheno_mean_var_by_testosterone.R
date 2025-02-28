library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)

# set up and load in files
pheno <- "height"
setwd("~/Research/GWAS-frontera/Phenotypes")
df_testosterone <- read.csv("pheno_testosterone.txt", sep="\t")
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t")
df_sex <- read.csv("sex_ids.txt", sep="\t")
colnames(df_pheno) <- c('FID','IID','pheno')

# merge dataframes and label sex
df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
df <- df[order(df$testosterone),]
df$sex[df$sex == 1] <- 'male'
df$sex[df$sex == 0] <- 'female'
df_m <- df[df$sex == 'male',]
df_f <- df[df$sex == 'female',]

# find intercept for overlaps
f_inter = mean(df_f$testosterone) + 2*sd(df_f$testosterone)
m_inter = mean(df_m$testosterone) - 2*sd(df_m$testosterone)

# histogram of testosterone frequency between sex
#ggplot(df, aes(x=testosterone, fill=sex)) +
#  geom_histogram(alpha=0.5, position='identity', bins=50) +
#  geom_vline(xintercept = f_inter, color='red') +
#  geom_vline(xintercept = m_inter, color='blue')

# create dataframe with the overlaps
overlap <- df_m[df_m$testosterone <= m_inter,]
overlap <- rbind(overlap, df_f[df_f$testosterone >= f_inter,])
nrow(overlap)
# remove overlapping from the non-overlapping df
df_m <- df_m[! df_m$IID %in% overlap$IID,]
df_f <- df_f[! df_f$IID %in% overlap$IID,]


### have a bin that have overlapping testosterone
### larger bins perhaps - have around 5-10 (other than that overlapping bin))

# mean and variance from each bin
bin_fun <- function(data, n, sex) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  for (i in 1:n) {
    bin <- data[cuts == levels(cuts)[i],]
    pheno_mean <- mean(bin$pheno)
    pheno_var <- var(bin$pheno)
    T_mean <- mean(bin$testosterone)
    results <- rbind(results, data.frame(Testosterone=T_mean, P_mean=pheno_mean, P_var=pheno_var, sex=sex))
  }
  return(results)
}

# resulting dataframes for graph
m_results <- bin_fun(df_m,10,'male')
f_results <- bin_fun(df_f,10,'female')
overlap_results_m <- bin_fun(overlap[overlap$sex=='male',],1,'male')
overlap_results_f <- bin_fun(overlap[overlap$sex=='female',],1,'female')
overlap_results <- rbind(overlap_results_m, overlap_results_f)
overlap_long <- melt(overlap_results, id.vars=c("Testosterone","sex"))
results <- rbind(m_results, f_results)
results_long <- melt(results, id.vars=c("Testosterone","sex"))
head(overlap_long)

# trendlines
m_over <- overlap_results$Testosterone[1]
f_over <- overlap_results$Testosterone[2]
x_max <- max(results$Testosterone)
mm <- lm(P_mean ~ Testosterone, data=results[results$sex == 'male',])
fm <- lm(P_mean ~ Testosterone, data=results[results$sex == 'female',])
mv <- lm(P_var ~ Testosterone, data=results[results$sex == 'male',])
fv <- lm(P_var ~ Testosterone, data=results[results$sex == 'female',])

trend_y <- function(m_model, f_model) {
  trendline <- NULL
  m_y1 <- m_model$coefficients[1] + (m_model$coefficients[2] * m_over)
  m_y2 <- m_model$coefficients[1] + (m_model$coefficients[2] * x_max)
  f_y1 <- f_model$coefficients[1] + (f_model$coefficients[2] * 0)
  f_y2 <- f_model$coefficients[1] + (f_model$coefficients[2] * f_over)
  trendline <- rbind(trendline,data.frame(x1=m_over, x2=x_max, y1=m_y1, y2=m_y2, sex='male'))
  trendline <- rbind(trendline,data.frame(x1=0, x2=f_over, y1=f_y1, y2=f_y2, sex='female'))
  return(trendline)
}

trend <- rbind(trend_y(mm,fm), trend_y(mv,fv))
trend$variable <- c('P_mean','P_mean','P_var','P_var')
trend
# plot
png( paste0(pheno,"_testosterone.png"), width=600, height=500, res=100)
labels <- c(P_mean="Phenotypic Mean", P_var="Phenotypic Variance")
ggplot(results_long, aes(x=Testosterone, y=value, color=sex)) +
  geom_point(size=2) +
  geom_point(data=overlap_long, aes(x=Testosterone, y=value), shape=1,stroke=1.5) +
  labs(title=paste0("Waist:Hip (BMI adj.)"," to Testosterone Levels"), x="Testosterone Level", y="Variance                                              Mean") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  facet_wrap(~variable, ncol=1, scales="free_y", labeller=labeller(variable=labels)) +
  geom_segment(data=trend, aes(x=x1, xend=x2,y=y1,yend=y2)) + 
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc = 0.7, label.y.npc=0.6, show.legend=FALSE) +
  theme_classic() + scale_color_npg()
dev.off()

