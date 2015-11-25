## This file generates random growth data for fitting

## See Pilling et al (2002) equation 3 for generating these data

## true average population parameters
library(plyr)
library(ggplot2)

sample.vbgf <- function(t.vec, Linf, k,  t0=0){
    lengths <- Linf*(1-exp(-k*(t.vec-t0)))
    lengths <- lengths*exp(rnorm(n=length(t.vec), mean=0, sd=sigma.obs)-sigma.obs^2/2)
    data.frame(ages=t.vec, lengths=lengths)
}
sample.ages <- function(n.ages) {sample(1:Ntime, size=n.ages, replace=FALSE)}
sample.lengths <- function(Nfish, n.ages){
    Linf.vec <- Linf.mean*exp(rnorm(n=Nfish, 0, sd=Linf.sigma)-Linf.sigma^2/2)
    k.vec <- k.mean*exp(rnorm(n=Nfish, mean=0, sd=k.sigma)-k.sigma^2/2)
    dat <- ldply(1:Nfish, function(i) cbind(fish=i, sample.vbgf(t.vec=sample.ages(n.ages), Linf=Linf.vec[i], k=k.vec[i])))
    saveRDS(dat, paste0('growth_data_',Nfish,'.RDS'))
    dat
}

dat.plot <- sample.lengths(Nfish=500, n.ages=5)
g <- ggplot(dat.plot, aes(ages, lengths))+geom_point(alpha=.5)
ggsave('plots/data_500.png', g, width=9, height=5)

dat <- sample.lengths(Nfish=50000, n.ages=5)

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
