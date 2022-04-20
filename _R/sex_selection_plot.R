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
load("fst_plot_whole_body_fat_mass.RData"); wbfm <- pointsf
load("fst_plot_waist_circ.RData"); waist <- pointsf
load("zscore_plot.RData")

# split results by ancestry
azj <- results[results$ANC == "Ashkenazi Jewish",]
fin <- results[results$ANC == "Finnish",]
nfe <- results[results$ANC == "Non-Finnish European",]

# calc ancestry zmean mean
head(resultsi)

### FST PLOT
# size 3x4
f1 <- ggplot(wbfm, aes(x=V, y=FST, weight=w)) +
  geom_point(color = "black", alpha= 0.1, size=1) +
  geom_smooth(method = "lm", se=F, color= "blue", size=0.5) +
  theme_classic() +
  scale_y_continuous(breaks=c(0,0.0015,0.003), labels=c("0","2e-3","3e-3")) +
  scale_x_continuous(breaks=c(0,0.003,0.006), labels=c("0","3e-3","6e-3")) +
  xlab("VGxSex (kg2)") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), axis.text=element_text(size=8))
f1


f2 <- ggplot(waist, aes(x=V, y=FST, weight=w)) +
  geom_point(color = "black", alpha= 0.1, size=1) +
  geom_smooth(method = "lm", se=F, color= "red", size=0.5) +
  theme_classic() +
  scale_y_continuous(breaks=c(0,0.00075, 0.0015), labels=c("0","7.5e-4","1.5e-3")) +
  scale_x_continuous(breaks=c(0,0.004,0.008), labels=c("0","4e-3","8e-3")) +
  xlab("VGxSex (cm2)") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), axis.text=element_text(size=8))
f2

fa <- ggplotGrob(f1) ; fb <- ggplotGrob(f2)
maxWidth = grid::unit.pmax(fa$widths[2:5], fb$widths[2:5])
fa$widths[2:5] <- as.list(maxWidth) ; fb$widths[2:5] <- as.list(maxWidth)
fplot <- grid.arrange(fa, fb, ncol=1, nrow=2)
annotate_figure(fplot, 
                left = text_grob("FST Between Males and Females", size=10, rot=90))

head(results)
### ZSCORE PLOT
# size 6x8
ggplot(results, aes(x=ZMEAN, y=TRAIT)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=CILOWER, xmax=CIUPPER), height=0.3, size=0.3) +
  geom_vline(xintercept = 0, linetype="longdash", color="#563f61", alpha=0.5, size=0.3) +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN), linetype="dashed", color="#b02e38", alpha=0.5, size=0.3) +
  coord_cartesian(xlim=c(-2,5)) +
  facet_wrap(~ANC,ncol=3) +
  theme_classic() +
  xlab("Z-score for Sexually-Antagonistic Selection") +
  theme(axis.title=element_text(size=12), axis.text=element_text(size=9),
        panel.grid.major.y = element_line(color="gray95", size=0.5))


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


