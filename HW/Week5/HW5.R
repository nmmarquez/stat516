rm(list=ls())

fnOpt <- function(par) {
    a <- par[1]
    b <- par[2]
    return( ( ((0.1 - qgamma(0.05,a,b))/0.1)^2 + ((10 - qgamma(0.95,a,b))/10)^2) )
}

opt <- optim(c(1,1),fnOpt,method="L-BFGS-B",lower=c(1e-100,1e-100),upper=c(1000,1000))


a <- opt$par[1]
b <- opt$par[2]
pgamma(10, shape=a, rate=b) - pgamma(.1, shape=a, rate=b)
pgamma(10, shape=a, rate=b) - pgamma(.1, shape=a, rate=b)
hist(rgamma(100000, a, b))
quantile(rgamma(100000, a, b), c(.05, .95))

Posterior <- rgamma(100000, 4+a, b+.25)
hist(Posterior)
quantile(Posterior, c(.025, .975))
qgamma(.025, 4+1, b+.25)
qgamma(.975, 4+1, b+.25)

snoqualmie <- readLines("./snoqualmie.txt")
snoqualmie <- gsub(" {2,}", " ", snoqualmie)
snoqualmie <- strsplit(snoqualmie, " ")
snoqualmie <- lapply(snoqualmie, function(years) as.numeric(years[-1]))

lapply(snoqualmie, length)

y <- scan("snoqualmie.txt")
nodays <- rep(c(365,365,365,366),9) # account for leap years
june <- matrix(0, nrow=36, ncol=30) # build empty data mtarix
daysum <- 31*3 + 30 + 29 # our start point is this number of days in the future
for (i in 1:36){
    if (i>1) daysum <- daysum + nodays[i-1] # if past the first year add 365|366
    june[i,] <- (y[daysum+1:30] > 0) * 1 # turn into indicators
}

n11 <- sum(june[,1:29] == 0 & june[,2:30] == 0)
n21 <- sum(june[,1:29] == 1 & june[,2:30] == 0)
n12 <- sum(june[,1:29] == 0 & june[,2:30] == 1)
n22 <- sum(june[,1:29] == 1 & june[,2:30] == 1)

# sanity check
sum(n11 + n12 + n21 + n22) == 36*29

p_hat_11 <- n11 / (n11 + n12)
p_hat_12 <- n12 / (n11 + n12)
p_hat_21 <- n21 / (n21 + n22)
p_hat_22 <- n22 / (n21 + n22)

P_hat <- rbind(c(p_hat_11, p_hat_12), c(p_hat_21, p_hat_22))
colnames(P_hat) <- c("dry", "wet")
row.names(P_hat) <- c("dry", "wet")

N_obs <- rbind(c(n11, n12), c(n21, n11))
colnames(N_obs) <- c("dry", "wet")
row.names(N_obs) <- c("dry", "wet")
P_stderr <- sqrt(P_hat * (1 - P_hat) / rowSums(N_obs))
P_hat - 1.96 * P_stderr
P_hat + 1.96 * P_stderr
P_hat

