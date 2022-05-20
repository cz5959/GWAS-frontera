#!/usr/bin/env Rscript

## DESKTOP R STUDIO ###

p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file and PCs
setwd("~/Research/GWAS-frontera/1000G")
phenotype <- read.table("EUR.height", header=T)
pcs <- read.table("EUR.eigenvec", header=F)

setwd("~/Research/GWAS-frontera/GWAS_Results/height/PRS")

# add the appropriate headers (6 PCs)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

# merge the files
pheno <- merge(phenotype, pcs, by=c("FID", "IID"))

# calculate the null model (model with PRS) using a linear regression 
# . refers to all variables 
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# R2 of the null model 
null.r2 <- summary(null.model)$r.squared

prs.result <- NULL
prs.residuals <- NULL
for(i in p.threshold){
    prs <- read.table(paste0("EUR_female_height_prs_2.",i,".profile"), header=T)
    # Merge prs with phenotype matrix - only taje FID, IID and PRS 
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # linear regression on Height with PRS and the covariates, ignoring the FID and IID
    model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
    prs.residuals <- rbind(prs.residuals, summary(model)$residuals)
}
prs_summary <- list("results" = prs.result, "residuals" = prs.residuals)
str(prs_summary)

# Best result is:
print(prs.result[which.max(prs.result$R2),])

#### PLOT ####
library(ggplot2)

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                    prs.result$print.p == 0] <-
    format(prs.result$P[!is.na(prs.result$print.p) &
                            prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
    # Specify that we want to print p-value on top of the bars
    geom_text(
        aes(label = paste(print.p)),
        vjust = -1.5,
        hjust = 0,
        angle = 45,
        cex = 4,
        parse = T
    )  +
    # Specify the range of the plot, *1.25 to provide enough space for the p-values
    scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
    # Specify the axis labels
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    # Draw a bar plot
    geom_bar(aes(fill = -log10(P)), stat = "identity") +
    # Specify the colors
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-4,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)
    ) +
    # Some beautification of the plot
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size =
                                        18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust =
                                    1)
    )
# save the plot
ggsave("EUR_female.height.bar.png", height = 7, width = 7)
