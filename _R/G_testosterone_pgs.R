library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)

# set up and load in files
pheno <- "SHBG"
title <- "SHBG"

# get PGS scores
setwd(paste0("~/Research/GWAS-frontera/GWAS_results/",pheno))
file_name <- list.files(pattern="both_sex_additive_")
df_both <- read.csv(file_name,sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))

# get phenotype 
setwd("~/Research/GWAS-frontera/Phenotypes")
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
df_sex <- read.csv("sex_ids.txt", sep="\t")

# testosterone pgs scores
setwd("~/Research/GWAS-frontera/GWAS_results/testosterone")
#df_testosterone <- read.csv("both_sex_additive_testosterone.1e-5.profile", 
#                            sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
## sex specific testosterone
df_T_f <- read.csv("female_additive_testosterone.1e-8.profile", 
                   sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_T_m <- read.csv("male_additive_testosterone.1e-5.profile", 
                   sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
df_testosterone <- merge(df_T_f,df_T_m, by="IID")

# merge dataframes - testosterone, sex, pheno, pgs scores f, pgs scores m
df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
## sex specific testosterone
df$SCORE.x <- ifelse(df$sex == 1, df$SCORE.x, df$SCORE.y)   
df <- df[-c(3)]

df <- merge(df,df_both,by='IID')
colnames(df) <- c("IID", "testosterone", "sex", "pheno", "SCORE")

# order by testosterone
df <- df[order(df$testosterone),]
# label then split by sex
df$sex[df$sex == 1] <- 'male'
df$sex[df$sex == 0] <- 'female'
df_m <- df[df$sex == 'male',]
df_f <- df[df$sex == 'female',]

## get n for bins
nrow(df_m)/10
nrow(df_f)/10

# get lm results from each testosterone bin
bin_fun <- function(data, n, sex) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  for (i in 1:n) {
    # linear regression
    bin <- data[cuts == levels(cuts)[i],]
    model <- lm(paste0("pheno ~ SCORE"), data = bin)
    beta <- model$coefficients[2]
    stderror <- summary(model)$coefficients[2,2]
    T_mean <- mean(bin$testosterone)
    results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Error=stderror, Sex=sex))
  }
  return(results)
}

# call function for overlap and nonoverlaps
m_results <- bin_fun(df_m,10,'male')
f_results <- bin_fun(df_f,10,'female')
results <- rbind(m_results, f_results)

# trendlines x range
x_min_m <- min(m_results$Testosterone); x_max_m <- max(m_results$Testosterone)
x_min_f <- min(f_results$Testosterone); x_max_f <- max(f_results$Testosterone)

# trendline boundary points
trend_y <- function(m_model, f_model) {
  trendline <- NULL
  m_inter <- m_model$coefficients[1]; m_slope <- m_model$coefficients[2]
  f_inter <- f_model$coefficients[1]; f_slope <- f_model$coefficients[2]
  m_y1 <- m_inter + (m_slope * x_min_m)
  m_y2 <- m_inter + (m_slope * x_max_m)
  f_y1 <- f_inter + (f_slope * x_min_f)
  f_y2 <- f_inter + (f_slope * x_max_f)
  mp <- summary(m_model)$coefficients[2,4]
  fp <- summary(f_model)$coefficients[2,4]
  trendline <- rbind(trendline,data.frame(x1=x_min_m, x2=x_max_m, y1=m_y1, y2=m_y2, p=mp, Sex='male'))
  trendline <- rbind(trendline,data.frame(x1=x_min_f, x2=x_max_f, y1=f_y1, y2=f_y2, p=fp, Sex='female'))
  return(trendline)
}

# linear regression for BETA
trend <- NULL
male_model <- lm(Beta ~ Testosterone, data=m_results)
female_model <- lm(Beta ~ Testosterone, data=f_results)
trend <- rbind(trend, trend_y(male_model, female_model))
text_size = 2.5
sex_label <- data.frame(label = c("female", "male"), Sex=c("female", "male"))
sex_label_y =8000
setwd("~/Research/GWAS-frontera/Supplementary Figures/Tpgs")

#png(file=paste0(pheno,"_Tpgs_sexspecific.png"), width=3.5, height=3, units="in", res=200)
ggplot(results, aes(x=Testosterone, y=Beta, color=Sex)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  labs(title=title, x="Testosterone PGS", y="Effect of PGS on Phenotype") +
  geom_segment(data=trend, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  theme_classic() + 
  theme(axis.text = element_text(size=8.5), axis.title = element_text(size=9), 
        plot.title=element_text(size=11), legend.position= "none", plot.margin = margin(20,20,20,20),
        strip.background = element_blank(), strip.text = element_blank()) +
  scale_color_manual(values=c("#d67629","#207335")) +
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc=0.0, label.y.npc = c(1,1), vjust=c(1,0.7),
           show.legend=FALSE, size=3) +
  geom_text(data=sex_label, aes(x=Inf, y=Inf, label=label), hjust=1.2, vjust=1.2, size=3.2) +
  facet_wrap(~Sex, ncol=1, scales="free")


#dev.off()

