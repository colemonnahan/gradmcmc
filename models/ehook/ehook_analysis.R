## This file reads in and processes the raw MCMC output, and then runs some
## analysis to validate the models and compare efficiency

par.names <- c("chain", "LP", "beta", "gamma","logcpue_mean", "logcpue_sd",
               "sigma_obs_mean", 'sigma_obs_sd', paste0('logcpue', 1:14),
                paste0('logsigma_obs', 1:14))


## Read in the data by algorithm and look at the indepent runs to make sure
## eveyrthing is working right
results.jags.ind <- readRDS(file='results/results.jags.ind.RDS')
xx <- results.jags.ind$BUGSoutput$sims.array
results.jags.ind <- do.call(rbind, lapply(1:dim(xx)[2], function(x) data.frame(chain=x, xx[,x,])))
names(results.jags.ind) <-
    gsub('.', '', x=names(results.jags.ind), fixed=TRUE)
chain <- results.jags.ind$chain
results.jags.ind$LP <- results.jags.ind$deviance*-2
results.jags.ind$deviance <- NULL
results.jags.ind <- results.jags.ind[, par.names]

## Stan
results.stan.ind <- readRDS(file='results/results.stan.ind.RDS')
temp <- get_sampler_params(results.stan.ind)
results.stan.tuning <- do.call(rbind, lapply(1:4, function(x)
    data.frame(chain=x, iteration=1:nrow(temp[[x]]), temp[[x]])))
stan.eps <- mean(subset(results.stan.tuning, iteration==max(results.stan.tuning$iteration))$stepsize__)
xx <- extract(results.stan.ind, permuted=FALSE)
results.stan.ind <- do.call(rbind, lapply(1:dim(xx)[2], function(x)
    data.frame(chain=x, xx[,x,])))
names(results.stan.ind) <-
    gsub('.', '', x=names(results.stan.ind), fixed=TRUE)
results.stan.ind$LP <- results.stan.ind$lp__
results.stan.ind$lp__ <- NULL
results.stan.ind <- results.stan.ind[, par.names]

## TMB results. these need to be processed manaully for the bounded
## parameters internally.
results.tmb.ind <- readRDS(file='results/results.tmb.ind.RDS')
results.tmb.ind$logcpue_mean=boundp(results.tmb.ind$logcpue_mean2, -5, 5)
results.tmb.ind$logcpue_sd=boundp(results.tmb.ind$logcpue_sd2, 0, 5)
results.tmb.ind$sigma_obs_mean=boundp(results.tmb.ind$sigma_obs_mean2, -5, 5)
results.tmb.ind$sigma_obs_sd=boundp(results.tmb.ind$sigma_obs_sd2, 0,5)
results.tmb.ind$beta= boundp(results.tmb.ind$beta2, 0,1)
results.tmb.ind$gamma=boundp(results.tmb.ind$gamma2, 0,1)
results.tmb.ind <- results.tmb.ind[, -grep('2', x=names(results.tmb.ind))]
names(results.tmb.ind) <-
    gsub('.', '', x=names(results.tmb.ind), fixed=TRUE)
names(results.tmb.ind)[2:15] <- paste0('logcpue', 1:14)
names(results.tmb.ind)[16:29] <- paste0('logsigma_obs', 1:14)
results.tmb.ind <- results.tmb.ind[, par.names]

## Make some plots to verify the models are the same
png('plots/ehook_qqplot_tmb_stan.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=c(1,1,1,1), oma=c(1,1,1,0), cex=.5, mgp=c(3,.5,0))
for(i in 1:ncol(results.stan.ind)){
    qqplot(results.tmb.ind[,i], results.stan.ind[,i], pch=16, col=rgb(1,0,0,.5))
    abline(a=0,b=1, lwd=2, col=rgb(0,0,0,.5))
    mtext(names(results.tmb.ind)[i], line=-2, cex=.5)
    mtext(names(results.stan.ind)[i], line=-1, cex=.5)
}
dev.off()
## Make some plots to verify the models are the same
png('plots/ehook_qqplot_tmb_jags.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=c(1,1,1,1), oma=c(1,1,1,0), cex=.5, mgp=c(3,.5,0))
for(i in 1:ncol(results.jags.ind)){
    qqplot(results.tmb.ind[,i], results.jags.ind[,i], pch=16, col=rgb(1,0,0,.5))
    abline(a=0,b=1, lwd=2, col=rgb(0,0,0,.5))
    mtext(names(results.tmb.ind)[i], line=-2, cex=.5)
    mtext(names(results.jags.ind)[i], line=-1, cex=.5)
}
dev.off()
## Make some plots to verify the models are the same
png('plots/ehook_qqplot_jags_stan.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=c(1,1,1,1), oma=c(1,1,1,0), cex=.5, mgp=c(3,.5,0))
for(i in 1:ncol(results.stan.ind)){
    qqplot(results.jags.ind[,i], results.stan.ind[,i], pch=16, col=rgb(1,0,0,.5))
    abline(a=0,b=1, lwd=2, col=rgb(0,0,0,.5))
    mtext(names(results.jags.ind)[i], line=-2, cex=.5)
    mtext(names(results.stan.ind)[i], line=-1, cex=.5)
}
dev.off()

## make set of plots for each model to verify
png('plots/ehook_tmb_acf.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.tmb.ind)) {
    acf(results.tmb.ind[,i], axes=FALSE);box()
    title(names(results.tmb.ind)[i], line=-1)
}
dev.off()
png('plots/ehook_stan_acf.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.stan.ind)) {
    acf(results.stan.ind[,i], axes=FALSE);box()
    title(names(results.stan.ind)[i], line=-1)
}
dev.off()
png('plots/ehook_jags_acf.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.jags.ind)) {
    acf(results.jags.ind[,i], axes=FALSE);box()
    title(names(results.jags.ind)[i], line=-1)
}
dev.off()


png('plots/ehook_tmb_trace.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.tmb.ind)){
    plot(results.tmb.ind[results.tmb.ind$chain==1,i], type='n', axes=FALSE,
         ylim=range(results.tmb.ind[,i])); box()
    for(j in unique(results.tmb.ind$chain)){
        lines(results.tmb.ind[results.tmb.ind$chain==j ,i ], col=j)
    }
    title(names(results.tmb.ind)[i], line=-1)
}
dev.off()
png('plots/ehook_stan_trace.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.stan.ind)){
    plot(results.stan.ind[results.stan.ind$chain==1,i], type='n', axes=FALSE,
         ylim=range(results.stan.ind[,i])); box()
    for(j in unique(results.stan.ind$chain)){
        lines(results.stan.ind[results.stan.ind$chain==j ,i ], col=j)
    }
    title(names(results.stan.ind)[i], line=-1)
}
dev.off()
png('plots/ehook_jags_trace.png', width=9, height=6, units='in', res=500)
par(mfrow=c(6,6), mar=.1*c(1,1,1,1))
for(i in 1:ncol(results.jags.ind)){
    plot(results.jags.ind[results.jags.ind$chain==1,i], type='n', axes=FALSE,
         ylim=range(results.jags.ind[,i])); box()
    for(j in unique(results.jags.ind$chain)){
        lines(results.jags.ind[results.jags.ind$chain==j ,i ], col=j)
    }
    title(names(results.jags.ind)[i], line=-1)
}
dev.off()

effectiveSize(results.tmb.ind)

ggplot(results.stan.tuning, aes(iteration, stepsize__, group=chain,
                                color=factor(chain))) + geom_line()
ggplot(results.stan.tuning, aes(iteration, n_leapfrog__, group=chain,
                                color=chain)) + geom_line()




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
