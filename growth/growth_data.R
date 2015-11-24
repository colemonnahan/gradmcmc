## This file generates random growth data for fitting

## See Pilling et al (2002) equation 3 for generating these data

## true average population parameters
library(plyr)
library(ggplot2)
## Global parameters
Linf.mean <- 50
k.mean <- .1
t0 <- 0
Linf.sigma <- .01
k.sigma <- .01
sigma.obs <- .1
Ntime <- 50

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

dat <- readRDS('growth_data_50000.RDS')

library(rstan)
Nfish <- 50
dat2 <- dat[dat$fish<=Nfish,]
stan.data <- list(Nfish=Nfish, Nobs=nrow(dat2), lengths=dat2$lengths,
                  fish=dat2$fish, ages=dat2$ages)
stan.init <- list(Linf_mean=Linf.mean, Linf_sigma=Linf.sigma,
                  k_mean=k.mean, k_sigma=k.sigma, sigma_obs=sigma.obs,
                  Linf=rep(Linf.mean, len=Nfish),
                  k=rep(k.mean, len=Nfish))
test <- stan(file='growth.stan', data=stan.data, init=list(stan.init),
             iter=500, warmup=100, chains=1)
data.frame(extract(test, permuted=TRUE)[,1,]

#' @param Linf An individual Linf, randomly generated
#' @param k An individual k, randomly generated
#' @param increment a single time increment, usually 1
#' @param Lt The length in the previous time period
vbgf.inc <- function(Lt, increment, Linf, k){
    Linf*(1-exp(-k*increment))+
        Lt*(exp(-k*increment))+ rnorm(1, 0, sigma.obs)
}

y <- vector(length=length(t.vec))
y[1] <- 0
for(tt in 2:(Ntime+ 1)){
    y[tt] <- vbgf.inc(y[tt-1], 1, Linf.mean, k.mean)
}
points(t.vec, y)
