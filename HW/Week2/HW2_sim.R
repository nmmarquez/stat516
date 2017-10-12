rm(list=ls())
set.seed(123)

N <- 1000000
p <- .75

# 2 independent binomial variables with p = .5
Y <- rbinom(N, 1, .5)
Z <- rbinom(N, 1, .5)

# make X a function of Y and Z susch that if Y == Z take 1 probability (p) and
# if Y != Z take 1 - p
eq_or_neq <- (Y + Z) %% 2
X <- rbinom(N, 1, p**(1-(eq_or_neq)) * (1-p)**(eq_or_neq))

RVs <- data.frame(Y, Z, X)
head(RVs)

summary(RVs)
summary(subset(RVs, Y == 1))
summary(subset(RVs, Z == 1))
summary(subset(RVs, Z == Y))
summary(subset(RVs, Z != Y))

