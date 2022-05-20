require("dplyr")
require("tidyr")
require("reshape2")
require("matrixStats")
require("ggplot2")
require("ggsci")
require("ggpubr")
require("gridExtra")

# load RData files (from Matt)
setwd("~/Research/GWAS-frontera/selection")
load("fst_plot_testosterone.1e-05.RData"); testosterone <- pointsf
load("fst_plot_protein_total.1e-05.RData"); protein <- pointsf
load("zscore_plot.1e-05.RData")
load("zscore_plot.1e-08.all.Rdata")

# split results by ancestry
azj <- results[results$ANC == "Ashkenazi Jewish",]
fin <- results[results$ANC == "Finnish",]
nfe <- results[results$ANC == "Non-Finnish European",]

nfe[nfe$TRAIT == "Total protein",c(6,7)]
summary(t_model)
# lm line
t_model <- lm(FST~V, testosterone, weight=w)
t_B
t_B <- summary(t_model)$coefficients[2]; t_yi <- summary(t_model)$coefficients[1]
p_model <- lm(FST~V, protein, weight=w)
p_B <- summary(p_model)$coefficients[2]; p_yi <- summary(p_model)$coefficients[1]
### FST PLOT
# size 3x4
f1 <- ggplot(testosterone, aes(x=V, y=FST, weight=w, size=w)) +
  geom_point(color = "black", alpha= 0.2) +
  #geom_segment(x=0,y=t_yi, xend=1.5e-3, yend=1.5e-3*t_B+t_yi, size= 0.5, color="blue") +
  geom_abline(slope=t_B, intercept=t_yi, size=0.5, color="blue") +
  theme_classic() +
  scale_y_continuous(breaks=c(0,5e-5,1e-4), labels=c("0","5e-5","1e-4"), limits = c(0,1e-4)) +
  scale_x_continuous(breaks=c(0,5e-4,1e-3, 1.5e-3), labels=c("0","5e-4","1e-3", "1.5e-3"), limits=c(0,1.5e-3)) +
  xlab("VGxSex (nmol/L)") +
  theme(axis.title.x=element_text(size=10), axis.text=element_text(size=8), axis.title.y=element_blank(),
        legend.position="none") +
  scale_size(range = c(0.5,3)) 
f1
max(testosterone$FST)
f2 <- ggplot(protein, aes(x=V, y=FST, weight=w,size=w)) +
  geom_point(color = "black", alpha= 0.2) +
  geom_abline(slope=p_B, intercept=p_yi, size=0.5, color="red") +
  theme_classic() +
  scale_y_continuous(breaks=c(0,4e-5, 8e-5), labels=c("0","4e-5","8e-5"), limits = c(0,9e-5)) +
  scale_x_continuous(breaks=c(0,5e-4,1e-3), labels=c("0","5e-4","1e-3"), limits=c(0,1.1e-3)) +
  xlab("VGxSex (g/L)") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), axis.text=element_text(size=8),
        legend.position="none") +
  scale_size(range = c(0.5,3)) 
f2

fa <- ggplotGrob(f1) ; fb <- ggplotGrob(f2)
maxWidth = grid::unit.pmax(fa$widths[2:5], fb$widths[2:5])
fa$widths[2:5] <- as.list(maxWidth) ; fb$widths[2:5] <- as.list(maxWidth)
fplot <- grid.arrange(fa, fb, ncol=1, nrow=2)
annotate_figure(fplot, 
                left = text_grob("FST Between Males and Females", size=10, rot=90))


### ZSCORE PLOT
# size 6x8
ggplot(results, aes(x=ZMEAN, y=TRAIT)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=CILOWER, xmax=CIUPPER), height=0.3, size=0.3) +
  geom_vline(xintercept = 0, linetype="longdash", color="#563f61", alpha=0.5, size=0.3) +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN), linetype="dashed", color="#b02e38", alpha=0.5, size=0.3) +
  coord_cartesian(xlim=c(-2.6,5.35)) +
  facet_wrap(~ANC,ncol=3) +
  theme_classic() +
  xlab("Z-score for Sexually-Antagonistic Selection") +
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=12), axis.text=element_text(size=9),
        panel.grid.major.y = element_line(color="gray95", size=0.5))
min(results$ZMEAN)
max(results$ZMEAN)
################3


z1 <- ggplot(azj, aes(x=ZMEAN, y=TRAIT)) +
  geom_point() +
  geom_errorbarh(aes(xmin=ZMEAN-ZSE, xmax=ZMEAN+ZSE), height=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", show.legend=F, color="black") +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN[1]), linetype="dashed", show.legend=F, color="#b02e38") +
  theme_classic() +
  ggtitle("Ashkenazi Jewish") +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=12), plot.title=element_text(size=14),
        axis.text=element_text(size=9))
z1

z2 <- ggplot(fin, aes(x=ZMEAN, y=TRAIT)) +
  geom_point() +
  geom_errorbarh(aes(xmin=ZMEAN-ZSE, xmax=ZMEAN+ZSE), height=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", show.legend=F, color="black") +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN[2]), linetype="dashed", show.legend=F, color="#b02e38") +
  theme_classic() +
  ggtitle("Finnish") +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), plot.title=element_text(size=14),
        axis.text=element_text(size=9))

z3 <- ggplot(nfe, aes(x=ZMEAN, y=TRAIT)) +
  geom_point() +
  geom_errorbarh(aes(xmin=ZMEAN-ZSE, xmax=ZMEAN+ZSE), height=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", show.legend=F, color="black") +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN[3]), linetype="dashed", show.legend=F, color="#b02e38") +
  theme_classic() +
  ggtitle("Non-Finnish European") +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), plot.title=element_text(size=14),
        axis.text=element_text(size=9))



zplot <- grid.arrange(z1, z2, z3, ncol=3, nrow=1)


