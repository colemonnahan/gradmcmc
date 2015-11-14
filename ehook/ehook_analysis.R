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
results.tmb.ind <- readRDS(file='results/results.tmb.ind.RDS')
## Look at pairs MLE covariance for long run.
mle2 <- list()
mle2$std <- 1*sqrt(diag(covar.tmb))
mle2$cor <- covar.tmb/ (sqrt(diag(covar.tmb)) %o% sqrt(diag(covar.tmb)))
mle2$nopar <- 34
mle2$names <- names(model.tmb$par)
mle2$est <- model.tmb.opt$par
xx <- results.tmb.ind
xx$seed <- xx$nll <- NULL
png(plots.file('ehook_tmb_ind_pairs.png'), width=9, height=6, units='in', res=1000)
pairs_ss(posterior=xx, mle=mle2, diag='acf', col=rgb(0,0,0,.5), pch='.')
dev.off()
rm(xx)
boundp <- function(x, min, max){
    min + (max-min)/(1+exp(-x))
}
results.tmb.ind$logcpue_mean=boundp(results.tmb.ind$logcpue_mean2, -5, 5)
results.tmb.ind$logcpue_sd=boundp(results.tmb.ind$logcpue_sd2, 0, 5)
results.tmb.ind$sigma_obs_mean=boundp(results.tmb.ind$sigma_obs_mean2, -5, 5)
results.tmb.ind$sigma_obs_sd=boundp(results.tmb.ind$sigma_obs_sd2, 0,5)
results.tmb.ind$beta= boundp(results.tmb.ind$beta2, 0,10)
results.tmb.ind$gamma=boundp(results.tmb.ind$gamma2, 0,1)
results.tmb.ind <- results.tmb.ind[, -grep('2', x=names(results.tmb.ind))]
names(results.tmb.ind) <-
    gsub('.', '_', x=names(results.tmb.ind), fixed=TRUE)
results.tmb.ind <- results.tmb.ind[, order(names(results.tmb.ind))]
coda::effectiveSize(results.tmb.ind)

png(plots.file('ehook_ind_acf.png'), width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:34) {
    acf(results.tmb.ind[,i], axes=FALSE);box()
    title(names(results.tmb.ind)[i], line=-1)
}
dev.off()

png(plots.file('ehook_ind_trace.png'), width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.tmb.ind)){
    plot(results.tmb.ind[,i], type='l', axes=FALSE, col=rgb(0,0,0,.5)); box()
    title(names(results.tmb.ind)[i], line=-1)
}
dev.off()

par(mfrow=c(5,7), mar=c(1,1,1,1), oma=c(1,1,1,1))
for(i in 1:ncol(results.tmb.ind)){
    qqplot(results.tmb.ind[,i], results2.jags[,i])
    abline(a=0,b=1)
    mtext(names(results.tmb.ind)[i], line=-1)
    mtext(names(results2.jags)[i], line=-2)
}
