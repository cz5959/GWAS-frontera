library("plotrix")
library("ggrepel")
library("ggplot2")
library("ggbreak")
library("gridExtra")
library("grid")

("~/Research/GWAS-frontera/Phenotypes")
df_vars <-  read.table("pheno_variances_bysex.txt", sep="\t",head=TRUE)
setwd("~/Research/GWAS-frontera/GWAS_Results")
df_weights = read.table("sum_mash_weights.txt", sep="\t",head=TRUE)

# formate dataframes
df_vars <- df_vars[df_vars$pheno %in% df_weights$phenotype,]
df_vars <- df_vars[order(df_vars$pheno),]
df_weights <- df_weights[order(df_weights$phenotype),]
row.names(df_vars) <- NULL ; row.names(df_weights) <- NULL

# get ratio and difference
df_vars$ratio_mtof <- df_vars$m_var / df_vars$f_var
df_weights$diff_mf <- df_weights$sum_weight_m - df_weights$sum_weight_f

outliers <- c("testosterone", "arm_fatfree_mass_R", "arm_fatfree_mass_L")
out_vars <- subset(df_vars, (df_vars$pheno %in% outliers))
out_weights <- subset(df_weights, (df_weights$phenotype %in% outliers))
out_x <- out_vars$ratio_mtof
out_y <- out_weights$diff_mf
main_vars <- subset(df_vars, !(df_vars$pheno %in% outliers))
main_weights <- subset(df_weights, !(df_weights$phenotype %in% outliers))
main_x <- main_vars$ratio_mtof
main_y <- main_weights$diff_mf

par(bty="n")
gap.plot(x, y, gap.axis="x", gap=c(4,27))
axis.break(1,4,breakcol="white",style="gap")
axis.break(1, 4*(1+0.02), breakcol="black", style="slash")
axis.break(2, 60*(1+0.02), breakcol="black", style="slash")

b <- element_blank()
p1 <- ggplot() + geom_point(aes(main_x,main_y)) + xlim(0,28) +
  geom_text(aes(main_x,main_y,label=main_weights$phenotype)) + 
  theme_classic() 
p2 <- ggplot() + geom_point(aes(out_x,out_y)) + xlim(0,28) + 
  scale_x_break(c(4,27)) +
  geom_text(aes(out_x,out_y,label=out_weights$phenotype)) + 
  theme_classic() + 
  theme(axis.title.x = b, axis.text.x = b, axis.line.x = b, axis.ticks.x=b)
grid.newpage()
grid.draw(rbind(ggplotGrob(p2),ggplotGrob(p1)))

p3 <- grid.arrange(p2,p1,heights=c(2,6))

plot(3:10, main = "Axis break test")
axis.break()
axis.break(2, 2.9, style = "zigzag")
twogrp <- c(rnorm(10) + 4, rnorm(10) + 20)
gap.plot(twogrp,gap = c(8,16), xlab = "Index", ylab = "Group values",
         main = "Two separated groups with gap axis break",
         col = c(rep(2, 10), rep(3, 10)), ytics = c(3, 5, 18, 20))

