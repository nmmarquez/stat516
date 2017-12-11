rm(list=ls())
library(HMM)

log_sum <- function(loga, logb){
    loga_ <- ifelse(loga >= logb, loga, logb)
    logb_ <- ifelse(loga >= logb, logb, loga)
    loga_ + log1p(exp(logb_ - loga_))
}

# Function for generating the TP matrix, Emission matrix, and inital probs
gen_params <- function(){
    P <- matrix(c(.95, .1, .05, .9), 2, 2)
    E <- matrix(c(rep(1/6, 6), rep(.1, 5), .5), 2, 6, byrow=T)
    v <- c(.5, .5)
    return(list(P=P, E=E, v=v))
}

params <- gen_params()

# Build a Markov object to test our implementation vs a known functional one 
HMMobj <- initHMM(1:2, 1:6, params$v, params$P, params$E)

# Function for simulating undelying and obsrved values
simulate_run <- function(n){
    params <- gen_params()
    H <- rep(0, n)
    H[1] <- sample(1:2, 1, prob=params$v) # initiate the start which is 50/50
    for(i in 2:n){
        H[i] <- sample(1:2, 1, prob=params$P[H[i-1],])
    }
    O <- sapply(1:n, function(i) sample(1:6, 1, prob=params$E[H[i],]))
    return(list(O=O, H=H))
}

# algorithm implementations
forward_algorithm <- function(obs, log_, scale_=FALSE){
    n <- length(obs)
    params <- gen_params()
    prob_seq <- matrix(0, nrow=2, ncol=n)
    prob_seq[,1] <- params$v * params$E[,obs[1]]
    if(log_){
        lP <- log(params$P)
        prob_seq[,1] <- log(prob_seq[,1])
        for(i in 2:n){
            prob_seq[,i] <- log(params$E[,obs[i]]) + 
                apply(lP + prob_seq[,i-1], 2, function(x) log_sum(x[1], x[2]))
        }
        prob_seq <- cbind(log(params$v), prob_seq)
    }
    else{
        if(scale_){
            prob_seq[,1] <- prob_seq[,1] / sum(prob_seq[,1])
        }
        for(i in 2:n){
            prob_seq[,i] <- params$E[,obs[i]] * c(prob_seq[,i-1] %*% params$P)
            if(scale_){
                prob_seq[,i] <- prob_seq[,i] / sum(prob_seq[,i])
            }
        }
        prob_seq <- cbind(params$v, prob_seq)
    }
    return(prob_seq)
}

backward_algorithm <- function(obs, log_, scale_=FALSE){
    n <- length(obs)
    params <- gen_params()
    prob_seq <- matrix(0, nrow=2, ncol=n+1)
    prob_seq[,n+1] <- c(1, 1)
    if(log_){
        lP <- t(log(params$P))
        prob_seq[,n+1] <- log(prob_seq[,n+1])
        for(i in n:1){
            temp <- prob_seq[,i+1] + log(params$E[,obs[i]])
            prob_seq[,i] <- apply(lP + temp, 2, function(x) log_sum(x[1], x[2]))
        }
    }
    else{
        for(i in n:1){
            prob_seq[,i] <- 
                c((prob_seq[,i+1] * params$E[,obs[i]]) %*% t(params$P))
            if(scale_){
                prob_seq[,i] <- prob_seq[,i] / sum(prob_seq[,i])
            }
        }
    }
    return(prob_seq)
}

fb_algorithm <- function(obs, log_, scale_=FALSE){
    forward_probs <- forward_algorithm(obs, log_, scale_)
    backward_probs <- backward_algorithm(obs, log_, scale_)
    if(log_){
        n <- length(obs) + 1
        mult <- forward_probs + backward_probs
        colSum <- apply(mult, 2, function(x) log_sum(x[1], x[2]))
        results <- sapply(1:n, function(x) mult[,x] - colSum[x])
    }
    else{
        mult <- forward_probs * backward_probs
        results <- apply(mult, 2, function(x) x / sum(x))
    }
    return(results)
}

set.seed(124)
# First We ant to test to make sure our implementation returns the same result 
# as the known working implementation
sims1 <- simulate_run(10) # simulate values
(log_forward_prob_package <- forward(HMMobj, sims1$O)) # calc log probs using HMM
(forward_prob_package <- exp(log_forward_prob_package)) # exponentiate to get raw probs
(log_forward_prob_self <- forward_algorithm(sims1$O, log_=T))
(forward_prob_self <- forward_algorithm(sims1$O, log_=F)) # run our function
all.equal(c(forward_prob_package), c(forward_prob_self[,-1])) # test for similarity
all.equal(c(log_forward_prob_package), c(log_forward_prob_self[,-1]))

# looks like ourforward algorithm works

# Run the same process with the backwards algorithm
log_backward_prob_package <- backward(HMMobj, sims1$O)
backward_prob_package <- exp(log_backward_prob_package)
backward_prob_self <- backward_algorithm(sims1$O, log_=FALSE)
log_backward_prob_self <- backward_algorithm(sims1$O, log_=TRUE)
all.equal(c(log_backward_prob_package), c(log_backward_prob_self[,-1]))
all.equal(c(backward_prob_package), c(backward_prob_self[,-1]))

# Now check if the log, original, and scaled algorithms include the same results
# for the forward-backward algorithm
all.equal(c(fb_algorithm(sims1$O, log_=F)), 
          c(exp(fb_algorithm(sims1$O, log_=T))))
all.equal(c(fb_algorithm(sims1$O, log_=F)), 
          c(fb_algorithm(sims1$O, log_=F, scale_=T)))

# All three algorithms give the same results at the end

# now we want to keep running our model until we get nonsensical prob estimates
# on average how large of a model do we fail on for the basic algorithm?
fail_counts <- 100
fail_obs <- rep(0, fail_counts)

for(i in 1:fail_counts){
    trail_length <- 400
    success <- TRUE
    while(success){
        trail_length <- trail_length + 1
        fail <- any(is.na(log(fb_algorithm(simulate_run(trail_length)$O, F))))
        success <- !fail
    }
    fail_obs[i] <- trail_length
    print(paste0("Finished trail number ", i, "!!!"))
}

# histogram of size of first noticed failed observation
png("./HistogramFailure.png")
hist(fail_obs, nclass=20, main="First Failed Observation Length for 100 Trails",
     xlab="Length")
dev.off()

# now try in log space
fail_counts <- 100
fail_obs <- rep(0, fail_counts)
for(i in 1:fail_counts){
    trail_length <- 400
    success <- TRUE
    while(success){
        if(trail_length > 1000){
            break
        }
        trail_length <- trail_length + 1
        fail <- any(is.na(fb_algorithm(simulate_run(trail_length)$O, T)))
        success <- !fail
    }
    if(trail_length > 1000){
        print("Models are not failing to evaluate up to around length 1000!")
        break
    }
    fail_obs[i] <- trail_length
}

# doesnt look like its breaking anymore with the log transform lets try scaled
fail_counts <- 100
fail_obs <- rep(0, fail_counts)
for(i in 1:fail_counts){
    trail_length <- 400
    success <- TRUE
    while(success){
        if(trail_length > 1000){
            break
        }
        trail_length <- trail_length + 1
        fail <- any(is.na(log(fb_algorithm(simulate_run(trail_length)$O, F, T))))
        success <- !fail
    }
    if(trail_length > 1000){
        print("Models are not failing to evaluate up to around length 1000!")
        break
    }
    fail_obs[i] <- trail_length
}

# Failures are also averted with the scaled form