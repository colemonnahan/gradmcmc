## This assumes the working directory is in the same folder as this file.
## main.dir <- paste0(getwd(), '/')
## load libraries and such
library(coda)
library(TMB)
library(ggplot2)
library(plyr)
library(rstan)
library(R2jags)
library(snowfall)
library(reshape2)
main.dir <- 'C:/Users/Cole/gradmcmc/'
ggwidth <- 8
ggheight <- 5

## Functions
results.file <- function(file) paste0(main.dir,'results/', file)
plots.file <- function(file) paste0(main.dir,'plots/', file)
figures.file <- function(file) paste0(main.dir,'paper/figures/', file)
#' @param d Dimension of matrix
#' @param covar The covariance for off diagonals. Assumed the same for all
#' pairs of parameters.
#' @return A Hessian matrix (inverse of covariance matrix).
make.covar <- function(d, covar){
    x <- matrix(covar, nrow=d, ncol=d)
    diag(x) <- 1
    x
}
## TMB bounding functions, copied from  ADMB
boundpinv <- function(x, min, max){
    -log( (max-min)/(x-min) -1)
}
boundp <- function(x, min, max){
    min + (max-min)/(1+exp(-x))
}
## Get cumulative minimum ESS as a percentage
cum.minESS <- function(df, breaks=5){
    x <- floor(seq(1, dim(df)[1], len=breaks))
    x <- x[x>1]
    y <- lapply(1:length(x), function(i){
       df2 <- df[1:x[i],,, drop=FALSE]
       ess <- data.frame(monitor(df2, warmup=0, print=FALSE))$n_eff
      data.frame(iteration=x[i], pct.ess=100*min(ess)/dim(df2)[1])
 })
    do.call(rbind, y)
}

#' @param model.name Character string for model name
#' @param model.jags A compiled JAGS model. This will be updated in the
#' function.
#' @param model.stan A compiled Stan model. This will be updated in the
#' function.
#' @param seeds A vector of seeds to run across (sequentially) for each
#' algorithm.
#' @param Nout The number of resulting iterations after thinning and warmup
#' @param L A vector of integer values for HMC.
#' @param stan.burnin The length of burnin for Stan
#' @param jags.burnin The length of burnin for JAGS
#' @param n.thin The thinning rate.
#' @param sink Whether to sink console output to file "sink_progress.txt"
#' to cleanup console output. Defaults to TRUE. Makes it easier to see the
#' progress on the console.
#' @return A list of two lists. adapt.list is the adaptive results from
#' Stan, and perf.list is the performance metrics for each run.
run.chains <- function(model.name, seeds, Nout, L,
                       model.jags, data.jags, inits.jags, params.jags,
                       model.stan, data.stan, inits.stan, params.stan,
                       n.burnin, n.thin, sink=TRUE){
    perf.list <- list()
    adapt.list <- list()
    k <- 1
    for(seed in seeds){
        ## Now run single long chains without thinning and timing to get
        ## performance (minESS/time) for each of the methods
        ## Run a long one to ensure good tuning
        if(sink){
            sink(file='sink_progress.txt', append=TRUE, type='output')
            on.exit(sink())
        }
        message(paste("Starting run at", Sys.time()))
        message(paste('Starting seed',seed))
        set.seed(seed)
        temp <- jags(data=data.jags, parameters.to.save=params.jags, inits=inits.jags,
                           model.file=model.jags, n.chains=1, DIC=FALSE,
                           n.iter=n.burnin+10, n.burnin=n.burnin, n.thin=1)
        ## Rerun with those tunings to get efficiency
        time.jags <-
            as.vector(system.time(results.jags <- update(temp, n.iter=n.thin*Nout, n.thin=n.thin))[3])
        sims.jags <- results.jags$BUGSoutput$sims.array
        perf.jags <- data.frame(rstan::monitor(sims=sims.jags, warmup=0, print=FALSE, probs=.5))
        perf.list[[k]] <-
            data.frame(model=model.name, platform='jags', seed=seed,
                       Npar=dim(sims.jags)[3], time=time.jags,
                       minESS=min(perf.jags$n_eff),
                       medianESS=median(perf.jags$n_eff),
                       Nsims=dim(sims.jags)[1])
        k <- k+1
        rm(sims.jags, perf.jags)
        ## Use adaptation for eps and diagonal covariances, but remove those
        ## samples and time it took
        message('Starting stan.nuts model')
        results.stan.nuts <-
            stan(fit=model.stan, data=data.stan, iter=n.thin*Nout+n.burnin,
                 warmup=n.burnin, chains=1, thin=n.thin, algorithm='NUTS',
                 init=inits.stan, seed=seed,
                 control=list(adapt_engaged=TRUE))
        time.stan.nuts <- get_elapsed_time(results.stan.nuts)[2]
        sims.stan.nuts <- extract(results.stan.nuts, permuted=FALSE)
        perf.stan.nuts <- data.frame(rstan::monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
        adapt.list[[k]] <- data.frame(model=model.name, alg='NUTS', seed=seed,
                                      Npar=dim(sims.stan.nuts)[3]-1,
                                      Nsims=dim(sims.stan.nuts)[1],
                                      get_sampler_params(results.stan.nuts))
        adapt.list[[k]]$iteration <- 1:nrow(adapt.list[[k]])
        perf.list[[k]] <-
            data.frame(model=model.name, platform='stan.nuts', seed=seed,
                       Npar=dim(sims.stan.nuts)[3]-1, time=time.stan.nuts,
                       minESS=min(perf.stan.nuts$n_eff),
                       medianESS=median(perf.stan.nuts$n_eff),
                       Nsims=dim(sims.stan.nuts)[1])
        k <- k+1
        rm(results.stan.nuts, sims.stan.nuts, perf.stan.nuts)
        for(LL in L){
            message(paste0('Starting stan.hmc',LL,' model'))
            results.stan.hmc <-
                stan(fit=model.stan, data=data.stan, iter=n.thin*Nout+n.burnin,
                     warmup=n.burnin, chains=1, thin=n.thin, algorithm='HMC', seed=seed,
                     init=inits.stan, control=list(adapt_engaged=TRUE, int_time=LL))
            time.stan.hmc <- get_elapsed_time(results.stan.hmc)[2]
            sims.stan.hmc <- extract(results.stan.hmc, permuted=FALSE)
            adapt.list[[k]] <- data.frame(model=model.name, alg=paste0('HMC',LL), seed=seed,
                                          Npar=dim(sims.stan.hmc)[3]-1,
                                          Nsims=dim(sims.stan.hmc)[1],
                                          get_sampler_params(results.stan.hmc))
            adapt.list[[k]]$iteration <- 1:nrow(adapt.list[[k]])
            perf.stan.hmc <- data.frame(rstan::monitor(sims=sims.stan.hmc, warmup=0, print=FALSE))
            perf.list[[k]] <-
                data.frame(model=model.name, platform=paste0('stan.hmc',LL), seed=seed,
                           Npar=dim(sims.stan.hmc)[3]-1, time=time.stan.hmc,
                           minESS=min(perf.stan.hmc$n_eff),
                           medianESS=median(perf.stan.hmc$n_eff),
                           Nsims=dim(sims.stan.hmc)[1])
            k <- k+1
            rm(results.stan.hmc, sims.stan.hmc, perf.stan.hmc)
        }
    }
    return(invisible(list(adapt.list=adapt.list, perf.list=perf.list)))
}

plot.model.results <- function(perf.list, adapt.list){
    perf <- do.call(rbind, do.call(rbind, perf.list))
    model.name <- as.character(perf$model[1])
    perf$samples.per.time <- perf$minESS/perf$time
    perf$minESS <- 100*perf$minESS/perf$Npar
    perf.long <- melt(perf, c('model', 'platform', 'seed', 'Npar', 'Nsims'))
    perf.long <- ddply(perf.long, .(platform, Npar, variable), mutate,
                       mean.value=mean(value))
    g <- ggplot(perf.long, aes(Npar, log(value), group=platform, color=platform)) +
        geom_jitter(position=position_jitter(width=1.5, height=0), alpha=.5) +
            geom_line(data=perf.long, aes(Npar, log(mean.value))) +
                facet_wrap('variable') + ggtitle("Performance Comparison")
    ggsave(paste0('plots/',model.name, '_perf_time.png'), g, width=ggwidth, height=ggheight)
    adapt <- do.call(rbind.fill, do.call(rbind, adapt.list))
    adapt$log.stepsize=log(adapt$stepsize__)
    adapt.long <- melt(adapt, c('model', 'alg', 'seed', 'Npar', 'Nsims', 'iteration'))
    g <- ggplot(adapt, aes(iteration, log.stepsize, group=seed, color=factor(seed))) +
        geom_line(alpha=.5, size=.5)+ facet_grid(Npar~alg, scales='fixed') +
            ggtitle('Step size adaptation')
    ggsave(paste0('plots/',model.name, '_adapt_eps.png'), g, width=ggwidth, height=ggheight)
    g <- ggplot(subset(adapt.long, alg=='NUTS' & variable %in% c("treedepth__", "n_divergent__")),
                aes(iteration, value, group=seed, color=factor(seed))) +
        geom_line(alpha=.5, size=.5)+ facet_grid(variable~Npar,
                            scales='free') + ggtitle('NUTS tuning metrics')
    ggsave(paste0('plots/',model.name, '_adapt_NUTS.png'), g, width=ggwidth, height=ggheight)
}
## ## Use this to make function to compare between models later
## ggplot(perf, aes(Npar, log(minESS), group=platform,
##                   color=platform))+geom_jitter(position=position_jitter(width=.5, height=0), alpha=.5) +
##                       geom_line(data=perf, aes(Npar, log(mean.minESS)))
## ggsave('plots/growth_perf_minESS.png', width=9, height=5)
## ggplot(perf, aes(Npar, log(medianESS), group=platform,
##                   color=platform))+geom_jitter(position=position_jitter(width=.5, height=0), alpha=.5) +
##                       geom_line(data=perf, aes(Npar, log(mean.medianESS)))
## ggsave('plots/growth_perf_medianESS.png', width=9, height=5)
## ggplot(perf, aes(Npar, log(samples.per.time), group=platform,
##                   color=platform))+geom_jitter(position=position_jitter(width=.5, height=0), alpha=.5) +
##                       geom_line(data=perf, aes(Npar, log(mean.samples.per.time)))
## ggsave('plots/growth_perf_time.png', width=9, height=5)



### ------------------------------------------------------------
### OLD FUNCTIONS
## ## wrappers to run chains and return ESS and other metrics
## run_hmc <- function(obj, nsim, eps, L, covar=NULL, seed=NULL, diag=FALSE){
##     if(!is.null(seed)) set.seed(seed)
##     x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="HMC", L=L,
##                     eps=eps, diagnostic=TRUE, covar=covar, Madapt=Madapt)
##     ## discard the warmup
##     par <- x1$par[-(1:Madapt),]
##     minESS <- min(as.vector(coda::effectiveSize(par)))
##     x2 <- data.frame(covar=!is.null(covar), tuning=eps, algorithm='hmc', L=L, seed=seed,
##                      time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
##                      perf=log10(x1$time/minESS))
##     if(diag) return(x1) else return(x2)
## }
## run_nuts <- function(obj, nsim, inits=NULL, covar=NULL, delta, seed=NULL,
##                      Madapt, diag=FALSE, max_doubling=4){
##     if(!is.null(seed)) set.seed(seed)
##     x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="NUTS", diagnostic=TRUE, max_doubling=max_doubling,
##                    covar=covar, delta=delta, Madapt=Madapt, params.init=inits)
##     ## discard the warmup
##     par <- x1$par[-(1:Madapt),]
##     minESS <- min(as.vector(coda::effectiveSize(par)))
##     x2 <- data.frame(covar=!is.null(covar), tuning=delta, algorithm='nuts',
##                      seed=seed, time=x1$time, minESS=minESS, acceptance=NA,
##                      perf=log10(x1$time/minESS))
##     if(diag) return(x1) else return(x2)
## }
## run_rwm <- function(obj, nsim, alpha, covar, seed=NULL, ...){
##     if(!is.null(seed)) set.seed(seed)
##     x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="RWM", diagnostic=TRUE,
##                    covar=covar, alpha=alpha, ...)
##     minESS <- min(as.vector(coda::effectiveSize(x1$par)))
##     x2 <- data.frame(covar=!is.null(covar), tuning=alpha, algorithm='rwm', seed=seed,
##                      time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
##                      perf=log10(x1$time/minESS))
##     return(x2)
## }
## make.trace <- function(df.thinned, model, string){
##     nrows <- ceiling(sqrt(ncol(df.thinned)))
##     png(paste0('plots/', model, '.trace.', string,'.png'), width=9, height=6,
##         units='in', res=500)
##     par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
##     for(i in 1:ncol(df.thinned)){
##         plot(df.thinned[,i], type='l', axes=FALSE,
##              ylim=range(df.thinned[,i]), col=rgb(0,0,0,.5)); box()
##         title(names(df.thinned)[i], line=-1)
##     }
##     dev.off()
## }
## make.acf <- function(df, model, string){
##     nrows <- ceiling(sqrt(ncol(df)))
##     png(paste0('plots/', model, '.acf.', string,'.png'), width=9, height=6,
##         units='in', res=500)
##     par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
##     for(i in 1:ncol(df)) {
##         acf(df[,i], axes=FALSE);box()
##         title(names(df)[i], line=-1)
##     }
##     dev.off()
## }
