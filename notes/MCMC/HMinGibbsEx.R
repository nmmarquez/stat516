rm(list=ls())
library(dplyr)

Y <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
E <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05,  1.05,  2.1, 10.5)
N <- length(Y)

plot(Y)
a0 <- 2 # hyperparamter alpha
b0 <- 1 # hyperparamter beta
cat("prior mean of lambda = ", a0/b0, "\n")
cat("prior var of lambda = ", a0/b0^2, "\n")
curve(dgamma(x, a0, b0), 0, 10, main="Prior Distribution of lambda", 
      ylab="Density", xlab="")

for(i in 1:N){
    a_i <- a0 + Y[i]
    b_i <- b0 + E[i]
    post_mean_i <- a_i/b_i
    post_var_i <- a_i/b_i^2
    title_ <- paste0("Posterior Density for y_", i)
    xlab_ <- paste0("Mean: ", round(post_mean_i, 4), 
                    ", Var: ", round(post_var_i, 4))
    curve(dgamma(x, a_i, b_i), 0, 4, main=title_, ylab="Density", xlab=xlab_)
}

# now use hyperprior to reduce parameter dependency 
ha <- .1 
hb <- 1
nchain <- 100000
burnin <- 10000
lambda.post <- matrix(1, nrow=N, ncol=nchain)
b0.post <- rep(.1, nchain)

for(i in 2:nchain){
    lambda.post[,i] <- lambda <- rgamma(N, a0 + Y, b0 + E)
    b0.post[i] <- b0 <- rgamma(1, N * a0 + ha, sum(lambda.post[,i]) + hb)
}

for(i in 1:N){
    title_ <- paste0("Posterior Density for y_", i)
    hist(lambda.post[i, burnin:nchain], main=title_, nclass=30)
}

hist(b0.post[burnin:nchain], main="Posterior of Beta", nclass=30)

# now use MH in gibbs samplling to put an exponential prior on alpha
