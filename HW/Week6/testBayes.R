rm(list=ls())
library(MASS)
library(tidyverse)

set.seed(123)
Y <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
E <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05,  1.05,  2.1, 10.5)
N <- length(Y)

nchain <- 10000
burnin <- 1000
lNames <- paste0("lambda ", 1:N)

proposalfunction <- function(param){
    return(rnorm(1, mean=param, sd=.2))
}

l0 <- 10
a0 <- 10
b0 <- 10
ha <- .1 # gamma alpha hyperprior for beta 
hb <- 1 # gamma beta hyperprior on beta

posteriorDF <- expand.grid(
    par = c(paste0("lambda ", 1:N), "a 1", "b 1"), 
    chain=1:nchain,
    val=NA) %>%
    as_tibble %>%
    spread("par", "val") %>%
    mutate_at(lNames, function(x) l0) %>%
    mutate(`a 1`=a0, `b 1`=b0)

for(i in 2:nchain){

    posteriorDF[i, lNames] <- lambda <- rgamma(N, a0 + Y, b0 + E)
    posteriorDF[i, "b 1"] <- b0 <- rgamma(1, N * a0 + ha, sum(lambda) + hb)
    
    astar <- proposalfunction(a0)
    # generate a probability of accepting that is g(p*)/g(pi)
    paccept <- prod(lambda)^(astar-a0) * b0^(N * (astar-a0)) *
        (gamma(astar)/gamma(a0))^-N * exp(-astar+a0)
    if (astar > 0 & runif(1) < paccept){
        a0 <- astar
    }

    posteriorDF[i, "a 1"] <- a0
}

posteriorDF %>%
    gather("par", "val", -chain) %>%
    mutate(idx=str_split(par, " ", simplify=TRUE)[,2]) %>%
    mutate(par=str_split(par, " ", simplify=TRUE)[,1]) %>%
    ggplot(aes(x=chain, y=val, color=idx, group=idx)) +
    geom_line(alpha=.7) +
    theme_classic() +
    facet_wrap(~par, scales="free_y", nrow = 2) +
    guides(color=FALSE)

posteriorDF %>%
    filter(chain > burnin) %>%
    ggplot(aes(`a 1`, `b 1`)) +
    geom_point() +
    theme_classic() +
    stat_density_2d(aes(color = stat(level))) +
    scale_color_distiller(palette = "Spectral") +
    guides(color=FALSE) +
    labs(x="Posterior of Alpha", y="Posterior of Beta") +
    ggtitle("Bivariate Posterior")

posteriorDF %>%
    filter(chain > burnin) %>%
    ggplot(aes(`a 1`, `b 1`)) +
    theme_classic() +
    stat_density_2d(aes(fill = stat(level)), geom="polygon") +
    scale_fill_distiller(palette = "Spectral") +
    labs(x="Posterior of Alpha", y="Posterior of Beta") +
    ggtitle("Bivariate Posterior")

