#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
wd <- paste("/scratch1/08005/cz5959/GWAS_Results/",args[1],sep="")
both_file <- paste("both_sex_all.",args[1],".glm.linear",sep="")
female_file <- paste("female_all.",args[1],".glm.linear",sep="")
male_file <- paste("male_all.",args[1],".glm.linear",sep="")

setwd(wd)
library("qqman", lib.loc="/work/08005/cz5959/frontera/R/x86_64-pc-linux-gnu-library/4.0/")

results_as <- read.table(both_file,sep="\t",head=FALSE,col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
jpeg( paste("manhattan_",args[1],"_both.jpeg",sep="") )
manhattan(results_as,chr="CHROM",bp="POS",p="P",snp="ID", main = paste("Manhattan Plot: ",toupper(args[1])," - Both Sex", sep=""), annotatePval=0.00000005 )

results_as <- read.table(female_file,sep="\t",head=FALSE,col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
jpeg( paste("manhattan_",args[1],"_female.jpeg",sep="") )
manhattan(results_as,chr="CHROM",bp="POS",p="P",snp="ID", main = paste("Manhattan Plot: ",toupper(args[1])," - Female", sep=""), annotatePval=0.00000005 )

results_as <- read.table(male_file,sep="\t",head=FALSE,col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"))
jpeg( paste("manhattan_",args[1],"_male.jpeg",sep="") )
manhattan(results_as,chr="CHROM",bp="POS",p="P",snp="ID", main = paste("Manhattan Plot: ",toupper(args[1])," - Male", sep=""), annotatePval=0.00000005 )

dev.off()  




