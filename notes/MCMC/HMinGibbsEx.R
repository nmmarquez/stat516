rm(list=ls())
library(dplyr, ggplot2)

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
    plot(density(lambda.post[i, burnin:nchain]), main=title_)
}

hist(b0.post[burnin:nchain], main="Posterior of Beta", nclass=30)

# now use MH in gibbs samplling to put an exponential prior on alpha
proposalfunction <- function(param){
    return(rnorm(1, mean=param, sd=.2))
}

lambda.post <- matrix(1, nrow=N, ncol=nchain)
b0.post <- rep(1, nchain)
a0.post <- rep(2, nchain)
a0 <- 6
b0 <- 1

for(i in 2:nchain){
    lambda.post[,i] <- lambda <- rgamma(N, a0 + Y, b0 + E)
    b0.post[i] <- b0 <- rgamma(1, N * a0 + ha, sum(lambda) + hb)
    
    astar <- proposalfunction(a0)
    # generate a probability of accepting that is g(p*)/g(pi)
    paccept <- prod(lambda)^(astar-a0) * b0^(N * (astar-a0)) *
        (gamma(astar)/gamma(a0))^-N * exp(-astar+a0)
    if (astar > 0 & runif(1) < paccept){
        a0 <- astar
    }
    a0.post[i] <- a0
}

plot(a0.post, type="l")
plot(b0.post, type="l")

for(i in 1:N){
    title_ <- paste0("Posterior Density for y_", i)
    plot(density(lambda.post[i, burnin:nchain]), main=title_)
}

hist(b0.post[burnin:nchain], main="Posterior of Beta", nclass=30)
hist(a0.post[burnin:nchain], main="Posterior of Alpha", nclass=30)
