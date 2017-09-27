rm(list=ls())
pacman::p_load(ggplot2)
set.seed(123)

# Question 1: Gamma-Poisson
M <- 1000
N <- 10000

alphas <- runif(M, 0.4, 2)
betas <- runif(M, 0.4, 2)

obs <- sapply(1:M, function(i) rgamma(N, shape=alphas[i], rate=betas[i]))
qplot(apply(obs, 2, mean) - (alphas / betas), xlab="Observed - Expected",
      main="Simulated Mean for Gamma Poisson")

qplot(apply(obs, 2, var) - (alphas / betas**2), xlab="Observed - Expected",
      main="Simulated Variance for Gamma Poisson")


# question 2 Bremaud 1.7.1
rm(list=ls())

Sum <- function(x){
    sapply(x, function(y) sum(0:y))
}

M <- 100000
maxN <- 20
Ns <- sample.int(maxN, M, replace=TRUE)
X <- t(sapply(1:M, function(i) sample.int(Ns[i], 2, replace=TRUE)))
Xmax <- apply(X, 1, max)

qplot(((Xmax ** 2 + Sum(Xmax-1)) / (2 * Xmax - 1)) - X[,1], 
      xlab="Observed - Expected", main="Simulated Conditional Means")
