#install.packages("mashr")
library(ashr)
library(mashr)

set.seed(1)
simdata = simple_sims(500,2,0.5)
data = mash_set_data(simdata$Bhat, simdata$Shat)

U.c = cov_canonical(data)  
  m = mash(data, Ulist= U.c, outputlevel = 1)

mp <- get_estimated_pi(m)

g <- get_fitted_g(m)

m2 = mash(data, g=get_fitted_g(m), fixg=TRUE)
m2 = mash(data, g=g, fixg=TRUE)
pm = get_pm(m2)


test <- function(x){
  y <- x+2
  y2 <- x+5
  return(list(y,y2))
}
results <- test(2)
results[[2]]
