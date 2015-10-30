## This file generates random growth data for fitting

## See Pilling et al (2002) equation 3 for generating these data

## true average population parameters
Linf.mean <- 50
k.mean <- .1
t0 <- 0
Nfish <- 10                             # number of individual fish
Ntime <- 50
Linf.sigma <- .01
k.sigma <- .01
obs.sigma <- 1
Linf.vec <- rnorm(n=Nfish, mean=Linf.mean, sd=Linf.sigma)
k.vec <- rnorm(n=Nfish, mean=k.mean, sd=k.sigma)
vbgf <- function(t.vec, Linf, k, t0=0) Linf*(1-exp(-k*(t.vec-t0)))
t.vec <- 0:Ntime
plot(t.vec, vbgf(t.vec, Linf.mean, k.mean), type='l')

#' @param Linf An individual Linf, randomly generated
#' @param k An individual k, randomly generated
#' @param increment a single time increment, usually 1
#' @param Lt The length in the previous time period
vbgf.inc <- function(Lt, increment, Linf, k){
    Linf*(1-exp(-k*increment))+
        Lt*(exp(-k*increment))+ rnorm(1, 0, obs.sigma)
}

y <- vector(length=length(t.vec))
y[1] <- 0
for(tt in 2:(Ntime+ 1)){
    y[tt] <- vbgf.inc(y[tt-1], 1, Linf.mean, k.mean)
}
points(t.vec, y)
