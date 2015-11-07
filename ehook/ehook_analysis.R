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

inits <- as.list(apply(results2.jags, 2, mean))
inits$LP <- NULL

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
