library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)
library(gridExtra)
library(grid)
library(dplyr)
library(boot)

# 1 male; 0 female
set.seed(1)

# phenotype names
setwd("~/Research/Phenotypes/")
pheno_names <- read.csv("pheno_names.txt", sep="\t")

# bootstrap lm R2
bootstrap_r2 <- function(df, v) {
  df_sample <- df[v,]
  model <- lm(pheno ~ SCORE, df_sample)
  bootstrap_R2 <- summary(model)$r.squared
  return(bootstrap_R2)
}

# function for retrieving lm values
get_lm_coefficients <- function(df, sex) {
  results <- NULL
  for (i in 0:1) {
    sex_df <- df[df$sex == i,]
    model <- lm(pheno ~ SCORE, sex_df); model <- summary(model)
    b <- model$coefficients[2]
    err <- model$coefficients[2,2]
    y_inter <- model$coefficients[1]
    boot_r2 <- boot(sex_df, bootstrap_r2, R=100)
    r2.se <- summary(boot_r2)$bootSE
    r2 <- model$r.squared
    
    results <- rbind(results, data.frame(Phenotype = pheno, Beta=b, Err=err, Y0=y_inter, R2=r2, R2.Err=r2.se, Sex=i))
  }
  return( results )
}

# get means from each bin
bin_fun <- function(data, n) {
  results <- NULL
  for (j in 0:1) {
    sex_df <- data[data$sex == j,]
    intervals = seq(0,nrow(sex_df),nrow(sex_df)/n)
    cuts <- cut(1:nrow(sex_df), breaks = intervals)
    for (i in 1:n) {
      # bin
      bin <- sex_df[cuts == levels(cuts)[i],]
      # get mean
      pheno_mean <- mean(bin$pheno)
      pheno_se <- sd(bin$pheno) / sqrt(nrow(bin))
      pgs_mean <- mean(bin$SCORE)
      results <- rbind(results, data.frame(Pheno=pheno_mean, Pheno_SE=pheno_se, Score=pgs_mean, Sex=j))
    }
  }
  return(results)
}

#slope_list <- NULL

for (i in 1:27) {
  pheno <- "hip_circ"
  pheno_title <- "Hip circumference"
  unit <- "cm"
#pheno <- pheno_names$Code[i]
#pheno_title <- pheno_names$Phenotype[i]
#unit <- pheno_names$Unit[i]
#print(pheno)

# load phenotype files
setwd("~/Research/Phenotypes")
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
df_sex <- read.csv("sex_ids.txt", sep="\t")
df_pheno <- merge(df_pheno, df_sex, by="IID")
colnames(df_pheno) <- c("IID", "pheno","sex")

# load pgs file
setwd(paste0("~/Research/GWAS-frontera/GWAS_results/",pheno))
file_name <- list.files(pattern="male_additive_")[2]
df_male <- read.csv(file_name, sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
file_name <- list.files(pattern="female_additive_")
df_female <- read.csv(file_name, sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))

# merge phenotype and pgs df
df_m <- merge(df_pheno, df_male, by="IID")
df_f <- merge(df_pheno, df_female, by="IID")

# order by score
df_m <- df_m[order(df_m$SCORE),]
df_f <- df_f[order(df_f$SCORE),]

# convert to score sd
df_m$SCORE <- df_m$SCORE / sd(df_m$SCORE)
df_f$SCORE <- df_f$SCORE / sd(df_f$SCORE)

# bin
m_results <- bin_fun(df_m, 10)
f_results <- bin_fun(df_f, 10)

# regress
m_lm <-  get_lm_coefficients(df_m, "male")
f_lm <-  get_lm_coefficients(df_f, "female")

# add to slope list
#slope_list <- rbind(rbind(slope_list, m_lm[-c(4)]), f_lm[-c(4)])

#plot
vdown = 1
label_m <- m_lm %>%
  mutate(Beta = sprintf("%.2f",Beta)) %>%
  mutate(Err = sprintf("%.2f",Err)) %>%
  mutate(R2 = sprintf("%.1f",R2*100)) %>%
  mutate(R2.Err = sprintf("%.1f",R2.Err*100))
label_f <- f_lm %>%
  mutate(Beta = sprintf("%.2f",Beta)) %>%
  mutate(Err = sprintf("%.2f",Err)) %>%
  mutate(R2 = sprintf("%.1f",R2*100)) %>%
  mutate(R2.Err = sprintf("%.1f",R2.Err*100))

setwd("~/Research/GWAS-frontera/Supplementary Figures/phenotype and pgs")
#png(file=paste0(pheno,"_pheno_pgs.png"), width=4, height=5, units="in", res=200)
#pdf(file=paste0(pheno,"_pheno_pgs.pdf"), width=3, height=4)
pdf(file=paste0(pheno,"_pheno_pgs2.pdf"), width=2, height=3)
p1 <- ggplot(data = m_results, aes(x=Score, y=Pheno, color=as.character(Sex))) +
  geom_point() +
  geom_abline(data=m_lm, aes(slope=Beta, intercept = Y0, color=as.character(Sex))) +
  theme_classic() +
  xlab("Male PGS SD") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), 
        axis.text = element_text(size=8.5), legend.position="none") + 
  scale_color_manual(values=c("#d67629", "#207335")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=0.05, vjust=vdown, size=3,color="#207335",
           label= bquote(slope == .(label_m$Beta[2])*"\u00b1"*.(label_m$Err[2])*", ")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=0.05, vjust=vdown+1, size=3,color="#207335",
           label=bquote( R^2 == .(label_m$R2[2])*"%\u00b1"*.(label_m$R2.Err[2])*"%" )) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=-1.1, vjust=vdown+9.8, size=3,color="#d67629",
           label= bquote(slope == .(label_m$Beta[1])*"\u00b1"*.(label_m$Err[1])*", ")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=-1.3, vjust=vdown+10.3, size=3,color="#d67629",
           label=bquote( R^2 == .(label_m$R2[1])*"%\u00b1"*.(label_m$R2.Err[1])*"%" )) 

p2 <- ggplot(data = f_results, aes(x=Score, y=Pheno, color=as.character(Sex))) +
  geom_point() +
  geom_abline(data=f_lm, aes(slope=Beta, intercept = Y0, color=as.character(Sex))) +
  theme_classic() +
  xlab("Female PGS SD") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), 
        axis.text = element_text(size=8.5), legend.position="none") + 
  scale_color_manual(values=c("#d67629", "#207335")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=0.05, vjust=vdown, size=3,color="#207335",
         label= bquote(slope == .(label_f$Beta[2])*"\u00b1"*.(label_f$Err[2])*", ")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=0.05, vjust=vdown+1, size=3,color="#207335",
           label=bquote( R^2 == .(label_f$R2[2])*"%\u00b1"*.(label_f$R2.Err[2])*"%" )) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=-1.1, vjust=vdown+9.8, size=3,color="#d67629",
           label= bquote(slope == .(label_f$Beta[1])*"\u00b1"*.(label_f$Err[1])*", ")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=-1.3, vjust=vdown+10.3, size=3,color="#d67629",
           label=bquote( R^2 == .(label_f$R2[1])*"%\u00b1"*.(label_f$R2.Err[1])*"%" )) 

p <- grid.arrange(p1, p2, ncol=1, nrow=2,
                  #top=textGrob(paste0("Sex-specific ", pheno_title, " and PGS"), gp=gpar(fontsize=12)),
                  left=textGrob(paste0(pheno_title, " (",unit,")"), gp=gpar(fontsize=10), rot=90))

print(p)
dev.off()

}

#pheno <- rep(pheno_names$Code,each=2)
#pheno_title <- rep(pheno_names$Phenotype,each=2)
#slope_list <- cbind(pheno, pheno_title, slope_list)
#setwd("~/Research/GWAS-frontera/PGS")
#write.table(slope_list, file = "sexspecific_pheno_pgs_lm.txt", sep="\t", row.names = F, quote = F)

# p <- ggplot(data = df, aes(x=Score, y=Pheno)) +
#   geom_point(aes(color=Sex)) +
#   geom_errorbar(aes(ymin=Pheno-Pheno_SE, ymax=Pheno+Pheno_SE, color=Sex), alpha= 0.5) +
#   geom_abline(slope=m_lm$Beta, intercept = m_lm$Y0, color="#207335") +
#   geom_abline(slope=f_lm$Beta, intercept = f_lm$Y0, color="#d67629") +
#   theme_classic() +
#   labs(title=paste0("Sex-specific ", pheno_title, " and PGS"),
#        x="PGS", y=pheno_title) +
#   theme(plot.title=element_text(size=12), axis.title=element_text(size=10), 
#         axis.text = element_text(size=8.5), legend.position = "none") + 
#   scale_color_manual(values = c("#d67629","#207335")) +
#   annotate("text", x=min(df$Score), y=max(df$Pheno), hjust=-0.1, vjust=1, size=3,
#            label=paste0("slope= ",round(m_lm$Beta,0), "\nR^2= ", round(m_lm$R2,2)), color="#207335") +
#   annotate("text", x=max(df$Score), y=min(df$Pheno), hjust=1, vjust=-0.5, size=3,
#            label=paste0("slope= ",round(f_lm$Beta,0), "\nR^2= ", round(f_lm$R2,2)), color="#d67629")
# print(p)
# ### OWN STUFF
# b_m <- slope_list[slope_list$Sex == "male",]
# b_f <- slope_list[slope_list$Sex == "female",]
# b_m$ratio_beta <- b_m$Beta / b_f$Beta
# 
# # get mash weights
# setwd("~/Research/GWAS-frontera/mash/")
# mash_weights <- read.csv("mash_weights.txt", sep="\t")
# # create column of difference of male and female mash weights
# mash_weights <- mash_weights %>%
#   mutate(diff_mf = sum_weight_m - sum_weight_f) %>%
#   select(c(1,6))
# 
# plot_df <- merge(b_m, mash_weights, by.x="pheno", by.y="phenotype")
# plot_df <- plot_df[c(1,2,7,8)]
# 
# test <- plot_df[plot_df$pheno != "creatinine",]
# ggplot(test, aes(diff_mf, ratio_beta)) +
#   geom_point() +
#   geom_label(aes(label=pheno_title))
  


