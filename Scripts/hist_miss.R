args <- commandArgs(trailingOnly = TRUE)
n <- args[1]
titleind <- sprintf("Histogram individual missingness Chr%s",n)
titlesnp <- sprintf("Histogram SNP missingness Chr%s",n)

indmiss<-read.table(file="plink2.smiss", header=TRUE)
snpmiss<-read.table(file="plink2.vmiss", header=TRUE)

# read data into R
# pdf("histsmiss_",n,".pdf") #indicates pdf format and gives title to file
pdf(sprintf("histsmiss_%s.pdf",n))
hist(indmiss[,5],main=titleind) #selects column 6, names header of file

pdf(sprintf("histvmiss_%s.pdf",n))
hist(snpmiss[,5],main=titlesnp)
dev.off() # shuts down the current device

