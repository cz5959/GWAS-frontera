#install.packages("mashr")
library(ashr)
library(mashr)

set.seed(1)
simdata = simple_sims(500,2,1)
data = mash_set_data(simdata$Bhat, simdata$Shat)
U.pca = cov_pca(data,5,subset=strong)

U.c = cov_canonical(data)
f_1 <- matrix(c(1,2,2,4),2,2)
f_2 <- matrix(c(1,3,3,9),2,2)
m_1 <- matrix(c(4,2,2,1),2,2)
m_2 <- matrix(c(9,3,3,1),2,2)
f_neg <- matrix(c(1,-2,-2,4),2,2)
m_neg <- matrix(c(4,-2,-2,1),2,2)
equal_opp <- matrix(c(1,-1,-1,1),2,2)
U.c[['f_1']] <- f_1
U.c[['f_2']] <- f_2
U.c[['m_1']] <- m_1
U.c[['m_2']] <- m_2
U.c[['f_neg']] <- f_neg
U.c[['m_neg']] <- m_neg
U.c[['equal_opp']] <- equal_opp

names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")
names(U.c)
names(U.ed) <- c("ED.PCA.1", "ED.PCA.2", "ED.tPCA")
m.c = mash(data, U.c)




x <- c(4, 1, 4, 5, 1, 6)
