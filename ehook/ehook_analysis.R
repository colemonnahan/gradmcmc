## This file reads in and processes the raw MCMC output, and then runs some
## analysis to validate the models and compare efficiency

## Read in the data by algorithm
results1.jags.list <- readRDS(file='results/results1.jags.list.RDS')
results1.jags <- as.data.frame(results1.jags.list$BUGSoutput$sims.array[,1,])
results1.jags$LP <- results1.jags$deviance*-2
results1.jags$deviance <- NULL
results1.jags <- results1.jags[, order(names(results1.jags))]
results2.jags.list <- readRDS(file='results/results2.jags.list.RDS')
results2.jags <- as.data.frame(results2.jags.list$BUGSoutput$sims.array[,1,])
results2.jags$LP <- results2.jags$deviance*-2
results2.jags$deviance <- NULL
results2.jags <- results2.jags[, order(names(results2.jags))]
effectiveSize(results2.jags)

## inits <- as.list(round(apply(results2.jags, 2, mean),3))
## inits$LP <- NULL
## dput(inits)

results1.stan.list <- readRDS(file='results/results1.stan.list.RDS')
results1.stan <- as.data.frame(extract(results1.stan.list, permuted=FALSE)[,1,])
results1.stan$LP <- results1.stan$lp__
results1.stan$lp__ <- NULL
results1.stan <- results1.stan[, order(names(results1.stan))]
effectiveSize(results1.stan)

sort(names(results1.stan))
sort(names(results1.jags))


par(mfrow=c(5,7), mar=c(1,1,1,1), oma=c(1,1,1,1))
for(i in 1:ncol(results1.stan)){
    qqplot(results1.stan[,i], results1.jags[,i])
    abline(a=0,b=1)
    mtext(names(results1.stan)[i], line=-1)
    mtext(names(results1.jags)[i], line=-2)
}


## TMB results. these need to be processed manaully for the bounded
## parameters internally.
results.tmb.independent <- readRDS(file='results/results.tmb.independent.RDS')
## remove the warmup samples
results.tmb.independent <- results.tmb.independent[-(1:500),]
boundp <- function(x, min, max){
    min + (max-min)/(1+exp(-x))
}
results.tmb.independent$logcpue_mean=boundp(results.tmb.independent$logcpue_mean2, -5, 5)
results.tmb.independent$logcpue_sd=boundp(results.tmb.independent$logcpue_sd2, 0, 5)
results.tmb.independent$sigma_obs_mean=boundp(results.tmb.independent$sigma_obs_mean2, -5, 5)
results.tmb.independent$sigma_obs_sd=boundp(results.tmb.independent$sigma_obs_sd2, 0,5)
results.tmb.independent$beta= boundp(results.tmb.independent$beta2, 0,10)
results.tmb.independent$gamma=boundp(results.tmb.independent$gamma2, 0,1)
results.tmb.independent <- results.tmb.independent[, -grep('2', x=names(results.tmb.independent))]
names(results.tmb.independent) <-
    gsub('.', '[', x=names(results.tmb.independent), fixed=TRUE)
results.tmb.independent <- results.tmb.independent[, order(names(results.tmb.independent))]
coda::effectiveSize(results.tmb.independent)



par(mfrow=c(6,6), mar=c(1,1,1,1))
for(i in 1:ncol(results.tmb.independent)){
   ## acf(results.tmb.independent[,i])
    plot(results.tmb.independent[,i], type='l')
}

par(mfrow=c(5,7), mar=c(1,1,1,1), oma=c(1,1,1,1))
for(i in 1:ncol(results.tmb.independent)){
    qqplot(results.tmb.independent[,i], results2.jags[,i])
    abline(a=0,b=1)
    mtext(names(results.tmb.independent)[i], line=-1)
    mtext(names(results2.jags)[i], line=-2)
}
