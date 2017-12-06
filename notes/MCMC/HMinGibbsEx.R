rm(list=ls())
library(dplyr)
library(ggplot2)
library(MASS)
library(RColorBrewer)

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
lambda.post1 <- matrix(1, nrow=N, ncol=nchain)
b0.post1 <- rep(.1, nchain)

for(i in 2:nchain){
    lambda.post1[,i] <- lambda <- rgamma(N, a0 + Y, b0 + E)
    b0.post1[i] <- b0 <- rgamma(1, N * a0 + ha, sum(lambda) + hb)
}

for(i in 1:N){
    title_ <- paste0("Posterior Density for y_", i)
    plot(density(lambda.post1[i, burnin:nchain]), main=title_)
}

hist(b0.post1[burnin:nchain], main="Posterior of Beta", nclass=30)

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

plot(density(b0.post[burnin:nchain]), main="Posterior of Beta")
plot(density(a0.post[burnin:nchain]), main="Posterior of Alpha")
hist(b0.post[burnin:nchain], nclass=30, main="Posterior of Beta", xlab="")
hist(a0.post[burnin:nchain], nclass=30, main="Posterior of Alpha", xlab="")

k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

z <- kde2d(a0.post[burnin:nchain], b0.post[burnin:nchain], n=50)

plot(a0.post[burnin:nchain], b0.post[burnin:nchain], 
     xlab=expression(alpha),  ylab=expression(beta), pch=19, cex=.4,
     main="Bivariate Posterior of Hyperparameters")
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)


# sample the joint posterior of alpha and beta
b0joint.post <- rep(1, nchain)
a0joint.post <- rep(2, nchain)
a0 <- 1
b0 <- 1

posterior <- function(a, b, theta=lambda.post[,N]){
    20 * a*log(b) - log(b^.9) + (-b*(1 + sum(theta)) - a) + log(gamma(a)^(-10)) + 
        log(prod(theta^(a-1)))
}

for(i in 2:nchain){
    astar <- proposalfunction(a0)
    bstar <- proposalfunction(b0)
    # generate a probability of accepting that is g(p*)/g(pi)
    paccept <- posterior(astar, bstar, theta=lambda.post[,i]) - 
        posterior(a0, b0, theta=lambda.post[,i])
    if (astar > 0 & bstar > 0 & log(runif(1)) < paccept){
        a0 <- astar
        b0 <- bstar
    }
    a0joint.post[i] <- a0
    b0joint.post[i] <- b0
}

plot(a0joint.post[burnin:nchain], type="l")
plot(b0joint.post[burnin:nchain], type="l")

hist(b0joint.post[burnin:nchain], nclass=30, 
     main="Posterior of Beta Given y", xlab="")
hist(a0joint.post[burnin:nchain], nclass=30, 
     main="Posterior of Alpha Given y", xlab="")

z <- kde2d(a0joint.post[burnin:nchain], b0joint.post[burnin:nchain], n=50)

plot(a0joint.post[burnin:nchain], b0joint.post[burnin:nchain], 
     xlab=expression(alpha),  ylab=expression(beta), pch=19, cex=.4,
     main="Bivariate Posterior of Hyperparameters given y")
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
