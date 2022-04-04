library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)

pheno_list <- c("albumin","arm_fatfree_mass_L","arm_fatfree_mass_R","bmi","calcium","creatinine","diastolicBP_auto",
                "eosinophil_perc","FVC_best","HbA1c","IGF1","lymphocyte_perc","protein_total","pulse_rate", "RBC_count",
                "systolicBP_auto","SHBG","urate","urea","waist_circ","waist_to_hip","weight","whole_body_fat_mass","wth_bmi_adj")
pheno_list <- c("waist_circ", "waist_to_hip", "weight", "whole_body_fat_mass", "pulse_rate","lymphocyte_perc", "albumin", "IGF1",
                "HbA1c", "FVC_best", "eosinophil_perc", "bmi")
pheno <- "hip_circ"
## male
for (pheno in pheno_list) {

  setwd(paste0("~/Research/GWAS-frontera/LDSC/",pheno))
  s_df <- read.csv(paste0(pheno,"_celltype_results.txt"), sep="\t", strip.white = TRUE,
                   colClasses = c("character", rep("numeric",2), rep("NULL",2)))
  colnames(s_df) <- c("Type","s_enrich","s_enrich_se")
  m_df <- read.csv(paste0(pheno,"_male_celltype_results.txt"), sep="\t", strip.white = TRUE,
                   colClasses = c("character", rep("numeric",2), rep("NULL",2)))
  colnames(m_df) <- c("Type","m_enrich","m_enrich_se")
  
  df <- merge(s_df, m_df, by="Type")
  max_e <- max(df$s_enrich+df$s_enrich_se,df$m_enrich+df$m_enrich_se)
  min_e <- min(df$s_enrich-df$s_enrich_se,df$m_enrich-df$m_enrich_se)
  
  df <- mutate(df, on = ifelse(abs(s_enrich-m_enrich) <= s_enrich_se | 
                                   abs(m_enrich-s_enrich) <= m_enrich_se, 1,2))
  setwd("~/Research/GWAS-frontera/LDSC/celltypes")
  png(file=paste0(pheno,"_male_sexdiff.png"), width=6,height=4,units="in",res=360)
  
  p <- ggplot(data=df, aes(x=m_enrich, y=s_enrich, color=as.factor(on))) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=s_enrich-s_enrich_se, ymax=s_enrich+s_enrich_se), width=0, show.legend = FALSE) +
    geom_errorbarh(aes(xmin=m_enrich-m_enrich_se, xmax=m_enrich+m_enrich_se), height=0, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="#2b62d9", size= 0.8) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(min_e,max_e)) + scale_y_continuous(limits = c(min_e,max_e)) +
    labs(title=paste0(pheno,": Sex-Diff v Male Celltype Partition"), x="Male Enrichment", y="M-F Enrichment") +
    geom_text_repel(aes(label=Type), max.overlaps = Inf, box.padding = 1.5) +
    scale_color_manual(values=c("#2b62d9","black"))
  print(p)
  
  dev.off()
}

## female
for (pheno in pheno_list) {
  setwd(paste0("~/Research/GWAS-frontera/LDSC/",pheno))
  s_df <- read.csv(paste0(pheno,"_sexdiff_celltype_results.txt"), sep="\t", strip.white = TRUE,
                   colClasses = c("character", rep("numeric",2), rep("NULL",2)))
  colnames(s_df) <- c("Type","s_enrich","s_enrich_se")
  f_df <- read.csv(paste0(pheno,"_female_celltype_results.txt"), sep="\t", strip.white = TRUE,
                   colClasses = c("character", rep("numeric",2), rep("NULL",2)))
  colnames(f_df) <- c("Type","f_enrich","f_enrich_se")
  
  df <- merge(s_df, f_df, by="Type")
  max_e <- max(df$s_enrich+df$s_enrich_se,df$f_enrich+df$f_enrich_se)
  min_e <- min(df$s_enrich-df$s_enrich_se,df$f_enrich-df$f_enrich_se)
  
  df <- mutate(df, on = ifelse(abs(s_enrich-f_enrich) <= s_enrich_se | 
                                 abs(f_enrich-s_enrich) <= f_enrich_se, 1,2))
  setwd("~/Research/GWAS-frontera/LDSC/celltypes")
  png(file=paste0(pheno,"_female_sexdiff.png"), width=6,height=4,units="in",res=360)
  
  p <- ggplot(data=df, aes(x=f_enrich, y=s_enrich, color=as.factor(on))) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=s_enrich-s_enrich_se, ymax=s_enrich+s_enrich_se), width=0, show.legend = FALSE) +
    geom_errorbarh(aes(xmin=f_enrich-f_enrich_se, xmax=f_enrich+f_enrich_se), height=0, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="#2b62d9", size= 0.8) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(min_e,max_e)) + scale_y_continuous(limits = c(min_e,max_e)) +
    labs(title=paste0(pheno,": Sex-Diff v Female Celltype Partition"), x="Female Enrichment", y="M-F Enrichment") +
    geom_text_repel(aes(label=Type), max.overlaps = Inf, box.padding =1) +
    scale_color_manual(values=c("#2b62d9","black"))
  print(p)
  
  dev.off()
}
