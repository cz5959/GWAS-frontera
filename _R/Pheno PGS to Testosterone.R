library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)

# set up and load in files
pheno <- "calcium"
title <- "Calcium"
setwd("~/Research/GWAS-frontera/Phenotypes")
df_testosterone <- read.csv("pheno_testosterone.txt", sep="\t", colClasses = c("NULL","integer","numeric"))
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
df_sex <- read.csv("sex_ids.txt", sep="\t")
colnames(df_pheno) <- c('IID','pheno')

# get PGS scores
setwd(paste0("~/Research/GWAS-frontera/GWAS_results/",pheno))

df_both <- read.csv("both_sex_additive_calcium.1e-5.profile",sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_add_f <- read.csv("female_additive_calcium.1e-5.profile",sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_add_m <- read.csv("male_additive_calcium.1e-5.profile",sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_mash_f <- read.csv("female_mash_calcium.1e-5.profile",sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_mash_m <- read.csv("male_mash_calcium.1e-5.profile",sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
score <- Reduce(function(x,y) merge(x,y, by='IID'), list(df_both, df_add_f, df_add_m, df_mash_f, df_mash_m))
colnames(score) <- c('IID', 'both', 'add_female', 'add_male', 'mash_female', 'mash_male')

# merge dataframes - testosterone, sex, pheno, pgs scores f, pgs scores m
df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
df <- merge(df,score,by='IID')

# order by testosterone
df <- df[order(df$testosterone),]
# label then split by sex
df$sex[df$sex == 1] <- 'male'
df$sex[df$sex == 0] <- 'female'
df_m <- df[df$sex == 'male',]
df_f <- df[df$sex == 'female',]

# find intercept for overlaps
f_inter = mean(df_f$testosterone) + 2*sd(df_f$testosterone)
m_inter = mean(df_m$testosterone) - 2*sd(df_m$testosterone)

# histogram of testosterone frequency between sex
#ggplot(df, aes(x=testosterone, fill=sex)) +
# geom_histogram(alpha=0.5, position='identity', bins=50) +
# geom_vline(xintercept = f_inter, color='red') +
# geom_vline(xintercept = m_inter, color='blue')

# create dataframe with the overlaps
overlap <- df_m[df_m$testosterone <= m_inter,]
overlap <- rbind(overlap, df_f[df_f$testosterone >= f_inter,])
# remove overlapping from the non-overlapping df
df_m <- df_m[! df_m$IID %in% overlap$IID,]
df_f <- df_f[! df_f$IID %in% overlap$IID,]

# get lm results from each testosterone bin
bin_fun <- function(data, n, sex) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  types <- c('both', 'add', 'mash')
  for (type in types) {
    if (type == 'both') {
      col_model <- 'both'
    } else {
      col_model <- paste0(type,"_",sex)
    }
    for (i in 1:n) {
      # linear regression
      bin <- data[cuts == levels(cuts)[i],]
      model <- lm(paste0("pheno ~ ",col_model), data = bin)
      beta <- model$coefficients[2]
      stderror <- summary(model)$coefficients[2,2]
      T_mean <- mean(bin$testosterone)
      # correlation
      corr <- cor.test(bin$pheno, unlist(bin[col_model]))
      corr_est <- corr$estimate
      corr_error <- corr$conf.int[2] - corr_est
      results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Error=stderror, Corr=corr_est, Corr_err=corr_error, Sex=sex, Type=type))
    }
  }
  return(results)
}

# call function for overlap and nonoverlaps
m_results <- bin_fun(df_m,10,'male')
f_results <- bin_fun(df_f,10,'female')
overlap_results_m <- bin_fun(overlap[overlap$sex=='male',],1,'male')
overlap_results_f <- bin_fun(overlap[overlap$sex=='female',],1,'female')
overlap_results <- rbind(overlap_results_m, overlap_results_f)
results <- rbind(m_results, f_results)

# trendlines x range
m_over <- overlap_results$Testosterone[1]
f_over <- overlap_results$Testosterone[4]
x_max <- max(results$Testosterone)

# trendline boundary points
trend_y <- function(m_model, f_model) {
  trendline <- NULL
  m_y1 <- m_model$coefficients[1] + (m_model$coefficients[2] * m_over)
  m_y2 <- m_model$coefficients[1] + (m_model$coefficients[2] * x_max)
  f_y1 <- f_model$coefficients[1] + (f_model$coefficients[2] * 0)
  f_y2 <- f_model$coefficients[1] + (f_model$coefficients[2] * f_over)
  mp <- summary(m_model)$coefficients[2,4]
  fp <- summary(f_model)$coefficients[2,4]
  trendline <- rbind(trendline,data.frame(x1=m_over, x2=x_max, y1=m_y1, y2=m_y2, p=mp, Sex='male'))
  trendline <- rbind(trendline,data.frame(x1=0, x2=f_over, y1=f_y1, y2=f_y2, p=fp, Sex='female'))
  return(trendline)
}

########################################################################

# linear regression for BETA
trend <- NULL
for (t in c('both', 'add', 'mash')) {
  result_sub <- results[results$Sex == 'male' & results$Type == t,]
  male_model <- lm(Beta ~ Testosterone, data=result_sub)
  result_sub <- results[results$Sex == 'female' & results$Type == t,]
  female_model <- lm(Beta ~ Testosterone, data=result_sub)
  trend <- rbind(trend, trend_y(male_model, female_model))
}
trend$Type <- c('both','both','add','add','mash','mash')

# plot phenotype to pgs linear regression
labels <- c(add="sex-specific additive", both="additive", mash="sex-specific covariance aware")
#png(paste0(pheno,"_pgs_testosterone.png"), width=700,height=700,res=100)
p <- ggplot(results, aes(x=Testosterone, y=Beta, color=Sex)) +
  geom_point(size=1.5) +
  geom_point(data=overlap_results, aes(x=Testosterone, y=Beta), shape=1, stroke=1, size=1.5) +
  geom_errorbar(aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  geom_errorbar(data= overlap_results, aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  labs(title=title, x="Testosterone Level", y="Phenotype Value to PGS Slope") +
  facet_wrap(~Type, ncol=1, scales="free_y", labeller=labeller(Type=labels)) +
  geom_segment(data=trend, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc=0.7, label.y.npc=0.90, show.legend=FALSE, size=3) +
  theme_classic() + 
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12), plot.title=element_text(size=16),
        legend.position= "none") +
  scale_color_manual(values=c("#d67629","#1d47a1"))

f_text <- data.frame(label = c("female", "", ""), Type = c("add", "both", "mash"))
m_text <- data.frame(label = c("male", "", ""), Type = c("add", "both", "mash"))

p + geom_text(data=f_text, x = 0.8, y = 11500, color="#d67629", aes(label=label), size=3.5, fontface=2) + 
  geom_text(data=m_text, x = 7, y = 11500, color="#1d47a1", aes(label=label), size=3.5, fontface=2)

#dev.off()

########## 
results_both <- results[results$Type == "both",]
trend_both <- trend[trend$Type == "both",]
overlap_results_both <- overlap_results[overlap_results$Type == "both",]
ggplot(results_both, aes(x=Testosterone, y=Beta, color=Sex)) +
  geom_point(size=1.5) +
  geom_point(data=overlap_results_both, aes(x=Testosterone, y=Beta), shape=1, stroke=1, size=1.5) +
  geom_errorbar(aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  geom_errorbar(data= overlap_results_both, aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  labs(title=title, x="Testosterone Level", y="Phenotype Value \nto PGS Slope") +
  geom_segment(data=trend_both, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc=0.6, label.y.npc=0.90, show.legend=FALSE, size=3) +
  theme_classic() + 
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12), plot.title=element_text(size=16),
        legend.position= "none", plot.margin = margin(20,20,20,20)) +
  scale_color_manual(values=c("#d67629","#207335")) +
  annotate("text", x=1, y=760, label="female", color="#d67629", size=3.5) +
  annotate("text", x=8, y=760, label="male", color="#207335", size=3.5)


########################################################################
# linear regression for CORR
trend <- NULL
for (t in c('both', 'add', 'mash')) {
  result_sub <- results[results$Sex == 'male' & results$Type == t,]
  male_model <- lm(Corr ~ Testosterone, data=result_sub)
  result_sub <- results[results$Sex == 'female' & results$Type == t,]
  female_model <- lm(Corr ~ Testosterone, data=result_sub)
  trend <- rbind(trend, trend_y(male_model, female_model))
}
trend$Type <- c('both','both','add','add','mash','mash')

# plot phenotype to pgs correlation
#png(paste0(pheno,"_pgs_testosterone_corr.png"), width=700,height=700,res=100)
ggplot(results, aes(x=Testosterone, y=Corr, color=Sex)) +
  geom_point(size=2) +
  geom_point(data=overlap_results, aes(x=Testosterone, y=Corr), shape=1, stroke=1.5) +
  geom_errorbar(aes(ymin=Corr-Corr_err, ymax=Corr+Corr_err), alpha= 0.4, show.legend = TRUE) +
  geom_errorbar(data= overlap_results, aes(ymin=Corr-Corr_err, ymax=Corr+Corr_err), alpha= 0.4) +
  labs(title=paste0(title," to PGS Correlation by Testosterone Levels"), x="Testosterone Level", y="Correlation") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  facet_wrap(~Type, ncol=1, scales="free_y", labeller=labeller(Type=labels)) +
  geom_segment(data=trend, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc=0.7, label.y.npc=0.15, show.legend=FALSE, size=3.5) +
  theme_classic() + scale_color_npg()
dev.off()

########################################################################
# weighted linear regression 
trend <- NULL
for (t in c('both', 'add', 'mash')) {7
  result_sub <- results[results$Sex == 'male' & results$Type == t,]
  male_model <- lm(Beta ~ Testosterone, data=result_sub, weight = 1 / result_sub$Error^2)
  result_sub <- results[results$Sex == 'female' & results$Type == t,]
  female_model <- lm(Beta ~ Testosterone, data=result_sub, weight = 1 / result_sub$Error^2)
  trend <- rbind(trend, trend_y(male_model, female_model))
}
trend$p <- paste0("p= ", round(trend$p,3))
trend$p[trend$p == 'p= 0'] <- "p < 0.001"
trend$Type <- c('both','both','add','add','mash','mash')

png(paste0(pheno,"_pgs_testosterone_weighted.png"), width=700,height=700,res=100)
ggplot(results, aes(x=Testosterone, y=Beta, color=Sex)) +
  geom_point(aes(size=Error), alpha=0.7) +
  geom_point(data=overlap_results, aes(x=Testosterone, y=Beta), shape=1, stroke=1.5) +
  geom_errorbar(aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4, show.legend = TRUE) +
  geom_errorbar(data= overlap_results, aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  labs(title=paste0(title," to PGS Weighted Regression \n by Testosterone Levels"), x="Testosterone Level", y="Phenotype to PGS") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), plot.title=element_text(size=20),
        legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  facet_wrap(~Type, ncol=1, scales="free_y", labeller=labeller(Type=labels)) +
  geom_segment(data=trend, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  geom_text(data=trend, aes(x=18, y=c(19000,21000,9000,10000,10000,11000), label=p), show.legend=FALSE) +
  theme_classic() + scale_color_npg()
dev.off()

########################


