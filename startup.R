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

#' Run Stan and JAGS models to compare efficiency.
#'
#' @param model Character string for model name
#' @param seeds A vector of seeds to run across (sequentially) for each
#' algorithm.
#' @param Nout The number of resulting iterations after thinning and warmup
#' @param Nthin The number to thin, defaults to 1.
#' @param lambda A vector of integrated time for HMC (eps times L)
#' @param delta A vector of target acceptance rates. Defaults to 0.8.
#' @param n.thin The thinning rate.
#' @param sink Whether to sink console output to trash file to cleanup
#' console output. Defaults to TRUE. Makes it easier to see the progress on
#' the console.
#' @param metric A vector of metrics to run across for HMC and NUTS. Must
#' be one of c("unit_e", "diag_e", "dense_e").
#' @return A list of two data frames. adapt is the adaptive results from
#' Stan, and perf is the performance metrics for each run.
run.chains <- function(model, seeds, Nout, Nthin=1, lambda, delta=.8,
                       metric='diag_e', data.jags, inits.jags,
                       params.jags, data.stan, inits.stan, sink.console=TRUE){
  model.jags <- paste0(model, '.jags')
  model.stan <- paste0(model, '.stan')
  if(Nthin != 1) stop("Nthin must be one, delta.final calculation breaks otherwise")
  Niter <- Nout*Nthin
  Nwarmup <- Niter/2
  ind.samples <- (Nwarmup/Nthin+1):Nout    # index of samples, excluding warmup
  ind.warmup <- 1:(Nwarmup/Nthin)             # index of warmup samples
  perf.list <- list()
  adapt.list <- list()

  ## Precompile Stan model so it isn't done repeatedly
  model.stan <-
    stan(file=model.stan, data=data.stan, iter=100,
         warmup=50, chains=1, thin=1, algorithm='NUTS',
         init=inits.stan, seed=1, verbose=FALSE,
         control=list(adapt_engaged=FALSE))
  k <- 1
  if(sink.console){
    sink(file='trash.txt', append=FALSE, type='output')
    on.exit(sink())
  }
  for(seed in seeds){
    message(paste('==== Starting seed',seed, 'at', Sys.time()))
    set.seed(seed)
    message('Starting JAGS model')
    time.jags.warmup <- as.vector(system.time(temp <-
      jags(data=data.jags, parameters.to.save=params.jags, inits=inits.jags,
           model.file=model.jags, n.chains=1, DIC=FALSE,
                      n.iter=Nwarmup+10, n.burnin=Nwarmup, n.thin=1)))[3]
    ## Rerun with those tunings to get efficiency
    time.jags.sampling <-
      as.vector(system.time(fit.jags <- update(temp, n.iter=Niter, n.thin=Nthin))[3])
    sims.jags <- fit.jags$BUGSoutput$sims.array
    perf.jags <- data.frame(rstan::monitor(sims=sims.jags, warmup=0, print=FALSE, probs=.5))
    perf.list[[k]] <-
      data.frame(platform='jags', seed=seed, Npar=dim(sims.jags)[3], delta.target=.5,
                 time.warmup=time.jags.warmup, time.sampling=time.jags.sampling,
                 minESS=min(perf.jags$n_eff), medianESS=median(perf.jags$n_eff),
                 Nsims=dim(sims.jags)[1])
    k <- k+1
    rm(sims.jags, perf.jags)
    ## Use adaptation for eps and diagonal covariances, but remove those
    ## samples and time it took
    message('Starting stan.nuts models')
    for(idelta in delta){
      for(imetric in metric){
        fit.stan.nuts <-
          stan(fit=model.stan, data=data.stan, iter=Niter+Nwarmup,
               warmup=Nwarmup, chains=1, thin=Nthin, algorithm='NUTS',
               init=inits.stan, seed=seed,
               control=list(adapt_engaged=TRUE, adapt_delta=idelta, metric=imetric))
        sims.stan.nuts <- extract(fit.stan.nuts, permuted=FALSE)
        perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
        adapt.nuts <- as.data.frame(get_sampler_params(fit.stan.nuts))
        adapt.list[[k]] <-
          data.frame(platform='NUTS', seed=seed,
                     Npar=dim(sims.stan.nuts)[3]-1,
                     Nsims=dim(sims.stan.nuts)[1],
                     delta.mean=mean(adapt.nuts$accept_stat__[ind.samples]),
                     delta.target=idelta,
                     eps.final=tail(adapt.nuts$stepsize__,1),
                     nsteps.mean=mean(adapt.nuts$n_leapfrog__[ind.samples]),
                     nsteps.median=median(adapt.nuts$n_leapfrog__[ind.samples]),
                     nsteps.sd=sd(adapt.nuts$n_leapfrog__[ind.samples]),
                     metric=imetric)
        perf.list[[k]] <-
          data.frame(platform='stan.nuts',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=tail(adapt.nuts$stepsize__,1),
                     Npar=dim(sims.stan.nuts)[3]-1,
                     time.sampling=get_elapsed_time(fit.stan.nuts)[2],
                     time.warmup=get_elapsed_time(fit.stan.nuts)[1],
                     minESS=min(perf.stan.nuts$n_eff),
                     medianESS=median(perf.stan.nuts$n_eff),
                     Nsims=dim(sims.stan.nuts)[1])
        k <- k+1
      }}
    rm(fit.stan.nuts, sims.stan.nuts, perf.stan.nuts, adapt.nuts)
    if(!is.null(lambda)){
      for(ilambda in lambda){
        message(paste0('Starting stan.hmc',ilambda,' models'))
        for(idelta in delta){
          for(imetric in metric){
            fit.stan.hmc <-
              stan(fit=model.stan, data=data.stan, iter=Niter+Nwarmup,
                   warmup=Nwarmup, chains=1, thin=Nthin, algorithm='HMC', seed=seed, init=inits.stan,
                   control=list(adapt_engaged=TRUE,
                     adapt_delta=idelta, metric=imetric, int_time=ilambda))
            sims.stan.hmc <- extract(fit.stan.hmc, permuted=FALSE)
            adapt.hmc <- as.data.frame(get_sampler_params(fit.stan.hmc))
            adapt.list[[k]] <-
              data.frame(platform=paste0('HMC',ilambda), seed=seed,
                         Npar=dim(sims.stan.hmc)[3]-1,
                         Nsims=dim(sims.stan.hmc)[1],
                         delta.target=idelta,
                         delta.mean=mean(adapt.hmc$accept_stat__[ind.samples]),
                         eps.final=tail(adapt.hmc$stepsize__,1),
                         metric=imetric, lambda=ilambda,
                         L=tail(adapt.hmc$stepsize__,1)*ilambda,
                         nsteps.mean=tail(adapt.hmc$stepsize__,1)*ilambda)
            perf.stan.hmc <- data.frame(rstan::monitor(sims=sims.stan.hmc, warmup=0, print=FALSE))
            perf.list[[k]] <-
              data.frame(platform=paste0('stan.hmc',ilambda),
                         seed=seed, delta.target=idelta, metric=imetric,
                         Npar=dim(sims.stan.hmc)[3]-1,
                         time.sampling=get_elapsed_time(fit.stan.hmc)[2],
                         time.warmup=get_elapsed_time(fit.stan.hmc)[1],
                         minESS=min(perf.stan.hmc$n_eff),
                         medianESS=median(perf.stan.hmc$n_eff),
                         Nsims=dim(sims.stan.hmc)[1])
            k <- k+1
            rm(fit.stan.hmc, sims.stan.hmc, perf.stan.hmc, adapt.hmc)
          }
        }
      }
    }
  }
  perf <- do.call(rbind.fill, perf.list)
  perf <- within(perf, {
                   time.total <- time.warmup+time.sampling
                   perf <- minESS/time.total
                 })
  adapt <- do.call(rbind.fill, adapt.list[!ldply(adapt.list, is.null)])
  perf$model <- adapt$model <- model
  return(invisible(list(adapt=adapt, perf=perf)))
}

#' Make plots comparing the performance of simulated data for a model.
## #'
## #' @details Makes plots with Npar (model complexity) on the x-axis and a
## #' variety of things on the y-axis.
## #' @param perf A data frame containing model performance data, as returned
## #' by run.chains
## #' @param adapt A data frame containing adaptation information from Stan,
## #' as returned by run.chains.
## #'
## #' @return Nothing. Makes plots in local folder of performance comparisons
## #' and adaptation results
## plot.simulated.results <- function(perf, adapt){
##     model.name <- as.character(perf$model[1])
##     perf$samples.per.time <- perf$minESS/perf$time
##     perf$minESS <- 100*perf$minESS/perf$Npar
##     perf.long <- melt(perf, c('model', 'platform', 'seed', 'Npar', 'Nsims'))
##     perf.long <- ddply(perf.long, .(platform, Npar, variable), mutate,
##                        mean.value=mean(value))
##     g <- ggplot(perf.long, aes(Npar, log(value), group=platform, color=platform)) +
##         geom_jitter(position=position_jitter(width=1.5, height=0), alpha=.5) +
##             geom_line(data=perf.long, aes(Npar, log(mean.value))) +
##                 facet_wrap('variable') + ggtitle("Performance Comparison")
##     ggsave(paste0('plots/',model.name, '_perf_time.png'), g, width=ggwidth, height=ggheight)
##     adapt$log.stepsize=log(adapt$stepsize__)
##     adapt.long <- melt(adapt, c('model', 'alg', 'seed', 'Npar', 'Nsims', 'iteration'))
##     g <- ggplot(adapt, aes(iteration, log.stepsize, group=seed, color=factor(seed))) +
##         geom_line(alpha=.5, size=.5)+ facet_grid(Npar~alg, scales='fixed') +
##             ggtitle('Step size adaptation')
##     ggsave(paste0('plots/',model.name, '_adapt_eps.png'), g, width=ggwidth, height=ggheight)
##     g <- ggplot(subset(adapt.long, alg=='NUTS' & variable %in% c("treedepth__", "n_divergent__")),
##                 aes(iteration, value, group=seed, color=factor(seed))) +
##         geom_line(alpha=.5, size=.5)+ facet_grid(variable~Npar,
##                             scales='free') + ggtitle('NUTS tuning metrics')
##     ggsave(paste0('plots/',model.name, '_adapt_NUTS.png'), g, width=ggwidth, height=ggheight)
## }

#' Make plots comparing the performance of empirical data for a model.
#'
#' @details Makes plots with delta (acceptance rate) on the x-axis and a
#' variety of things on the y-axis.
#' @param perf A data frame containing model performance data, as returned
#' by run.chains
#' @param adapt A data frame containing adaptation information from Stan,
#' as returned by run.chains.
#'
#' @return Nothing. Makes plots in local folder of performance comparisons
#' and adaptation results
plot.empirical.results <- function(perf, adapt){
    model.name <- as.character(perf$model[1])
    perf.long <-
      melt(perf, id.vars=c('platform', 'seed', 'delta.target', 'metric'),
           measure.vars=c('time.warmup', 'time.sampling', 'minESS',
             'perf'))
    perf.long$seed2 <- as.factor(with(perf.long, paste(seed, metric, sep="_")))
    g <- ggplot(perf.long, aes(delta.target, log(value), group=seed2, color=metric))+
      geom_line() + geom_point() + facet_grid(variable~platform, scales='free_y') + xlim(0,1)
    ggsave(paste0('plots/',model.name, '_perf.png'), g, width=ggwidth, height=ggheight)
    adapt.long <- melt(adapt, id.vars=c('platform', 'seed', 'delta.target', 'metric'),
                       measure.vars=c('eps.final', 'delta.mean', 'nsteps.mean'))
    adapt.long$seed2 <- as.factor(with(adapt.long, paste(seed, metric, sep="_")))
    g <- ggplot(adapt.long, aes(delta.target, value, group=seed2, color=metric))+
      geom_line() + facet_grid(variable~platform, scales='free_y') + xlim(0,1)
    ggsave(paste0('plots/',model.name, '_adapt.png'), g, width=ggwidth, height=ggheight)
  }


#' Verify the models are the same (coding errors) by plotting QQ plots of
#' each parameter between Stan and JAGS.
#'
#' @param sims.stan A data frame of MCMC samples from a single chain of
#' a Stan run, using extract(fit.stan).
#' @param sims.jags A data frame of MCMC samples from a single chain of a
#' JAGS run.
#' @return ggplot object invisibly. Also makes plot in folder 'plots' in
#' current working directory.
plot.model.comparisons <- function(sims.stan, sims.jags){
    ## Clean up names so they match exactly
    names(sims.stan) <- gsub('\\.', '', x=names(sims.stan))
    names(sims.jags) <- gsub('\\.', '', x=names(sims.jags))
    sims.stan$lp__ <- sims.jags$deviance <- NULL
    par.names <- names(sims.stan)
    sims.jags <- sims.jags[,par.names]
    ## Massage qqplot results into long format for ggplot
    qq <- ldply(par.names, function(i){
                    temp <- as.data.frame(qqplot(sims.jags[,i], sims.stan[,i], plot.it=FALSE))
                    return(cbind(i,temp))
                })
    g <- ggplot(qq, aes(x,y))+ geom_point(alpha=.5) +
        geom_abline(slope=1, col='red') + facet_wrap('i', scales='free') +
            xlab('jags')+ ylab('stan')
    ggsave('plots/model_comparison.png', width=9, height=5)
    return(invisible(g))
}
#' Run models with thinning to get independent samples which are then used
#' to verify the posteriors are the same, effectively checking for bugs
#' between models before doing performance comparisons
#'
verify.models <- function(params.jags, model.jags, inits.jags, data.jags,
                          params.stan, model.stan, inits.stan, data.stan,
                          Niter, Nthin){
    Nwarmup <- Niter/2
    fit.jags <- jags(data=data.jags, inits=inits.jags, param=params.jags,
                     model.file=model.jags, n.chains=1, n.burnin=Nwarmup, n.iter=Niter,
                     n.thin=Nthin)
    sims.jags <- data.frame(fit.jags$BUGSoutput$sims.matrix)
    array.jags <- fit.jags$BUGSoutput$sims.array
    fit.stan <- stan(file=model.stan, data=data.stan, iter=Niter, chains=1,
                     warmup=Nwarmup, thin=Nthin, init=inits.stan)
    sims.stan <- data.frame(extract(fit.stan, permuted=FALSE)[,1,])
    jags.ess <- data.frame(monitor(sims=fit.jags$BUGSoutput$sims.array, warmup=0, print=FALSE, probs=.5))$n_eff
    stan.ess <- data.frame(monitor(sims=extract(fit.stan, permuted=FALSE), warmup=0, print=FALSE, probs=.5))$n_eff
    plot.model.comparisons(sims.stan, sims.jags)
    print(rbind(jags.ess, stan.ess))
}
make.trace <- function(df.thinned, model, string){
    nrows <- ceiling(sqrt(ncol(df.thinned)))
    png(paste0('plots/', model, '.trace.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df.thinned)){
        plot(df.thinned[,i], type='l', axes=FALSE,
             ylim=range(df.thinned[,i]), col=rgb(0,0,0,.5)); box()
        title(names(df.thinned)[i], line=-1)
    }
    dev.off()
}
make.acf <- function(df, model, string){
    nrows <- ceiling(sqrt(ncol(df)))
    png(paste0('plots/', model, '.acf.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df)) {
        acf(df[,i], axes=FALSE);box()
        title(names(df)[i], line=-1)
    }
    dev.off()
}
