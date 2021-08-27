#install.packages("mashr")
library(ashr)
library(mashr)

set.seed(1)
simdata = simple_sims(500,2,0.5)
str(simdata)
head(simdata$Shat)
data = mash_set_data(simdata$Bhat, simdata$Shat)
  


