rm(list=ls())
library(HMM)

simulate_run <- function(n){
    P <- matrix(c(.95, .1, .05, .9), 2, 2)
    E <- matrix(c(rep(1/6, 6), rep(.1, 5), .5), 2, 6, byrow=T)
    v <- c(.5, .5)
    H <- rep(0, n)
    H[1] <- sample(1:2, 1, prob=v) # initiate the start which is 50/50
    for(i in 2:n){
        H[i] <- sample(1:2, 1, prob=P[H[i-1],])
    }
    O <- sapply(1:n, function(i) sample(1:6, 1, prob=E[H[i],]))
    return(list(O=O, H=H))
}

# algorithm implementations
forward_algorithm <- function(obs, E, P, v){
    n <- length(obs)
    prob_seq <- matrix(0, nrow=nrow(E), ncol=n)
    prob_seq[,1] <- v * E[,obs[1]]
    prob_seq[,1] <- prob_seq[,1] / sum(prob_seq[,1])
    for(i in 2:n){
        prob_seq[,i] <- E[,obs[i]] * c(prob_seq[,i-1] %*% P)
        prob_seq[,i] <- prob_seq[,i] / sum(prob_seq[,i])
    }
    prob_seq <- cbind(v, prob_seq)
    colnames(prob_seq) <- NULL
    return(prob_seq)
}

backward_algorithm <- function(obs, E, P, v, states){
    n <- length(obs)
    prob_seq <- matrix(0, nrow=nrow(E), ncol=n+1)
    prob_seq[,n+1] <- c(1, 1)
    for(i in n:1){
        prob_seq[,i] <- c((prob_seq[,i+1] * E[,obs[i]]) %*% t(P))
        prob_seq[,i] <- prob_seq[,i] / sum(prob_seq[,i])
    }
    return(prob_seq)
}

baumWelchSelf <- function(obs, niter=2000, states=1:2, emissions=1:6){
    # initialize parameters
    nStates <- length(states)
    nEmissions <- length(emissions)
    P <- matrix(1/nStates, nrow=nStates, ncol=nStates)
    diag(P) <- 1/nStates * (nStates:1)*2 
    P <- t(apply(P, 1, function(x) x / sum(x)))
    E <- matrix(1/nEmissions, nrow=nStates, ncol=nEmissions)
    v <- rep(1/nStates, nStates)
    n <- length(obs)
    
    for(q in 1:niter){
        a_it <- forward_algorithm(obs, E, P, v)[,2:(n+1)]
        b_it <- backward_algorithm(obs, E, P, v)[,2:(n+1)]
        e_itj <- array(0, dim=c(nStates, n-1, nStates))
        gamma_it <- a_it * b_it
        gamma_it <- apply(gamma_it, 2, function(x) x / sum(x))
        for(t in 1:(n-1)){
            for(i in 1:nStates){
                for(j in 1:nStates){
                    e_itj[i,t,j] <- a_it[i,t] * b_it[j,t+1] * P[i,j] * E[j,obs[t+1]]
                }
            }
            e_itj[,t,] <- e_itj[,t,] / sum(e_itj[,t,])
        }
        v <- gamma_it[,1]
        P <- t(apply(apply(e_itj, c(1,3), sum), 1, function(x) x/ sum(x)))
        for(i in 1:nStates){
            for(k in 1:nEmissions){
                E[i,k] <- sum((obs==k) * gamma_it[i,]) / sum(gamma_it[i,])
            }
        }
    print(paste0("Finished update number: ", q))
    }
    return(list(P=P, E=E, v=v, a=a_it, b=b_it, gamma=gamma_it, epsilon=e_itj))
}

set.seed(123)
sim1 <- simulate_run(5000)
est <- baumWelchSelf(sim1$O, 100)
HMMobj <- initHMM(
    1:2, 1:6, c(.5, .5), 
    matrix(c(.8, 1/3, .2, 2/3), 2, 2), 
    matrix(c(rep(1/6, 6), rep(1/6, 6)), 2, 6, byrow=T))

estPackage <- baumWelch(HMMobj, sim1$O)
estPackage$hmm$transProbs
est$P
estPackage$hmm$emissionProbs
est$E
abs(estPackage$hmm$emissionProbs - est$E) 
est$v

# read in file to run analysis on
weather <- c(t(as.matrix(read.table("./set3.dat", sep=""))))
estDat <- baumWelchSelf(weather, niter=100, emissions = 1:16)
HMMobj <- initHMM(
    1:2, 1:16, rep(1/2, 2),
    matrix(c(.8, 1/3, .2, 2/3), 2, 2), 
    matrix(c(rep(1/16, 16), rep(1/16, 16)), 2, 16, byrow=T))
estDatPackage <- baumWelch(HMMobj, weather)
estDatPackage$hmm$transProbs
estDat$P
estDatPackage$hmm$emissionProbs
estDat$E
abs(estDatPackage$hmm$emissionProbs - estDat$E)
estDat$E

estDat$a[,1000:1010]
estDat$b[,1:10]
tail(t(estDat$b), n=10)
