require("dplyr")
require("tidyr")
require("ashr")
require("mashr")

setwd("~/Research/GWAS-frontera/mash/simulation/nonnull")
load("mash_100_nonnull.RData")

# 64% female, 18% neither, 18% male
set.seed(1)
num <- 1:nrow(mash_BETA)
f <- sample(num, length(num)*0.64)
m <- sample(num[-f], length(num)*0.18)
b <- num[-c(f,m)]

mash_BETA$female[f] <- mash_BETA$female[f] * 2
mash_BETA$male[m] <- mash_BETA$male[m] * 2



mash_func <- function(mash_BETA, mash_SE) {
  ## mash
  mash_BETA <- as.matrix(mash_BETA) ; mash_SE <- as.matrix(mash_SE)
  data <- mash_set_data(mash_BETA, mash_SE, zero_Shat_reset = .Machine$double.eps)
  # set up covariance matrices
  U.c = cov_canonical(data)
  corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
  effect = c(1.5,2,3)
  for (c in corr) {
    for (e in effect) {
      U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
      U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
    }
  }
  U.c[['equal_-0.25_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
  U.c[['equal_-0.5_1']] <- matrix(c(1,-0.5,-0.5,1),2,2)
  U.c[['equal_-0.75_1']] <- matrix(c(1,-0.75,-0.75,1),2,2)
  U.c[['equal_-1_1']] <- matrix(c(1,-1,-1,1),2,2)
  names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")
  
  m.c <- mash(data,U.c)
  mixture <- get_estimated_pi(m.c)
  return(m.c)
}

m.c <- mash_func(mash_BETA, mash_SE)

write.table(mixture, file="mash_1000_64_18.txt",row.names = FALSE, sep="\t")



a <- matrix(1:12, ncol=3, nrow=4)
a <- c(1:8)
b <- c(1,3,5,7)
c <- c(2,4,6,8)
a[1:3]
a
b
c(colSums(a[,c(1,2)]), colSums(a[,c(1,2)]))




library("MASS")
n = 100
d = 2
mu <- c(0,0)
sigma_1 <- matrix(c(1,1,1,2), 2, 2)
Beta <- mvrnorm(10, mu=mu, Sigma=sigma_1)
Beta
Beta <- rbind(Beta, mvrnorm(30000*16, mu=mu, Sigma=sigma_2))

Beta[,1]
cov(Beta)

a[b] <- a[b]*2
a


