## This file generates random growth data for fitting

## See Pilling et al (2002) equation 3 for generating these data

## true average population parameters
library(plyr)
library(ggplot2)


library(rstan)





### ------------------------------------------------------------


## #' @param Linf An individual Linf, randomly generated
## #' @param k An individual k, randomly generated
## #' @param increment a single time increment, usually 1
## #' @param Lt The length in the previous time period
## vbgf.inc <- function(Lt, increment, Linf, k){
##     Linf*(1-exp(-k*increment))+
##         Lt*(exp(-k*increment))+ rnorm(1, 0, sigma.obs)
## }
## y <- vector(length=length(t.vec))
## y[1] <- 0
## for(tt in 2:(Ntime+ 1)){
##     y[tt] <- vbgf.inc(y[tt-1], 1, Linf.mean, k.mean)
## }
## points(t.vec, y)
