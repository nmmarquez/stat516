rm(list=ls())

likfunc <- function(x){
    b <- exp(x)
    (.9 - (-exp(-10*b) + exp(-.1*b)))**2
}

test <- optim(0, likfunc, method="Brent", lower=-20, upper=20)
b <- exp(test$par)
pgamma(10, shape=1, rate=b) - pgamma(.1, shape=1, rate=b)
pgamma(10, shape=1, rate=b) - pgamma(.1, shape=1, rate=b)

Posterior <- rgamma(100000, 4+1, b+.25)
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
    print(paste(daysum+1, daysum+30))
    june[i,] <- (y[daysum+1:30] > 0) * 1 # turn into indicators
}

n00 <- sum(june[,1:29] == 0 & june[,2:30] == 0)
n10 <- sum(june[,1:29] == 1 & june[,2:30] == 0)
n01 <- sum(june[,1:29] == 0 & june[,2:30] == 1)
n11 <- sum(june[,1:29] == 1 & june[,2:30] == 1)

# sanity check
sum(n00 + n01 + n10 + n11) == 36*29

p_hat_00 <- n00 / (n00 + n01)
p_hat_01 <- n01 / (n00 + n01)
p_hat_10 <- n10 / (n10 + n11)
p_hat_11 <- n11 / (n10 + n11)

P_hat <- rbind(c(p_hat_00, p_hat_01), c(p_hat_10, p_hat_11))
colnames(P_hat) <- c("dry", "wet")
row.names(P_hat) <- c("dry", "wet")

N_obs <- rbind(c(n00, n01), c(n10, n11))
colnames(N_obs) <- c("dry", "wet")
row.names(N_obs) <- c("dry", "wet")
P_stderr <- sqrt(P_hat * (1 - P_hat) / rowSums(N_obs))
P_hat - 1.96 * P_stderr
P_hat + 1.96 * P_stderr
P_hat
