rm(list=ls())
pacman::p_load(ggplot2)
set.seed(123)

# Question 1: Gamma-Poisson
M <- 1000
N <- 100000

alphas <- runif(M, 0.4, 2)
betas <- runif(M, 0.4, 2)

obs <- sapply(1:M, function(i) 
    rpois(N, rgamma(N, shape=alphas[i], rate=betas[i])))


qplot(apply(obs, 2, mean) - (alphas / betas), xlab="Observed - Expected",
      main="Simulated Mean for Gamma Poisson")

qplot(apply(obs, 2, var) - ((alphas * (betas + 1)) / betas**2), xlab="Observed - Expected",
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


# question 3 
rm(list=ls())

M <- 10000
N <- 10000

intensities <- sapply(1:M, function(x) runif(N, 0, 130))
rvariables <- sapply(1:M, function(x) rexp(N, intensities[,x]))
minvals <- apply(rvariables, 2, min)
minindex <- apply(rvariables, 2, function(x) which(x == min(x))) 

hist(minindex, main="", xlab="Index Value")
hist(minvals, main="", xlab="Min Value")

cor(minindex, minvals)

# lecture notes  
rm(list=ls())

m <- 1000
alleals <- sample(0:1, 2*m, replace=T)
years <- 5000

allealyears <- matrix(0, 2*m, years)
dim(allealyears)

allealyears[,1] <- alleals

for(y in 2:years){
    allealyears[,y] <- sample(allealyears[,y-1], 2*m, replace=T)
}

tail(colSums(allealyears))