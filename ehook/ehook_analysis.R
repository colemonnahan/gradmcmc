## This file reads in and processes the raw MCMC output, and then runs some
## analysis to validate the models and compare efficiency

## Read in the data by algorithm and look at the indepent runs to make sure
## eveyrthing is working right
results.jags.ind <- readRDS(file='results/results.jags.ind.RDS')
xx <- results.jags.ind$BUGSoutput$sims.array
results.jags.ind <- do.call(rbind, lapply(1:dim(xx)[2], function(x) data.frame(chain=x, xx[,x,])))
names(results.jags.ind) <-
    gsub('.', '', x=names(results.jags.ind), fixed=TRUE)
chain <- results.jags.ind$chain
results.jags.ind$LP <- results.jags.ind$deviance*-2
results.jags.ind$chain <- results.jags.ind$deviance <- NULL
results.jags.ind <- results.jags.ind[, order(names(results.jags.ind))]
results.jags.ind <- cbind(chain, results.jags.ind)

results.stan.ind <- readRDS(file='results/results.stan.ind.RDS')
xx <- extract(results.stan.ind, permuted=FALSE)
results.stan.ind <- do.call(rbind, lapply(1:dim(xx)[2], function(x)
    data.frame(chain=x, xx[,x,])))
names(results.stan.ind) <-
    gsub('.', '', x=names(results.stan.ind), fixed=TRUE)
chain <- results.stan.ind$chain
results.stan.ind$LP <- results.stan.ind$lp__
results.stan.ind$chain <- results.stan.ind$lp__ <- NULL
results.stan.ind <- results.stan.ind[, order(names(results.stan.ind))]
results.stan.ind <- cbind(chain, results.stan.ind)

results.tmb.ind <- readRDS(file='results/results.tmb.ind.RDS')
chain <- results.tmb.ind$chain
results.tmb.ind$chain <-  NULL
names(results.tmb.ind)
## TMB results. these need to be processed manaully for the bounded
## parameters internally.
results.tmb.ind$logcpue_mean=boundp(results.tmb.ind$logcpue_mean2, -5, 5)
results.tmb.ind$logcpue_sd=boundp(results.tmb.ind$logcpue_sd2, 0, 5)
results.tmb.ind$sigma_obs_mean=boundp(results.tmb.ind$sigma_obs_mean2, -5, 5)
results.tmb.ind$sigma_obs_sd=boundp(results.tmb.ind$sigma_obs_sd2, 0,5)
results.tmb.ind$beta= boundp(results.tmb.ind$beta2, 0,10)
results.tmb.ind$gamma=boundp(results.tmb.ind$gamma2, 0,1)
results.tmb.ind <- results.tmb.ind[, -grep('2', x=names(results.tmb.ind))]
names(results.tmb.ind) <-
    gsub('.', '', x=names(results.tmb.ind), fixed=TRUE)
results.tmb.ind <- results.tmb.ind[, order(names(results.tmb.ind))]
results.tmb.ind <- cbind(chain, results.tmb.ind)

par(mfrow=c(6,6), mar=c(1,1,1,1), oma=c(1,1,1,1))
for(i in 1:ncol(results.stan.ind)){
    qqplot(results.stan.ind[,i], results.jags.ind[,i], pch='.', col='red')
    abline(a=0,b=1, lwd=2)
    mtext(names(results.stan.ind)[i], line=-1)
    mtext(names(results.jags.ind)[i], line=-2)
}

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
    plot(results.tmb.ind[results.tmb.ind$chain==1,i], type='n', axes=FALSE, col=rgb(0,0,0,.5));
    box()
    for(j in unique(results.tmb.ind$chain)){
        lines(results.tmb.ind[results.tmb.ind$chain==j ,i ], col=j)
    }
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
