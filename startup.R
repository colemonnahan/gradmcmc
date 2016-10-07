## Load libraries, functions, and global variables
library(coda)
library(ggplot2)
library(plyr)
library(rstan)
library(R2jags)
library(reshape2)
library(MASS)                           # has mvrnorm()
ggwidth <- 8
ggheight <- 5
results.file <- function(file) paste0(main.dir,'results/', file)

#' Run Stan and JAGS models to compare efficiency.
#'
#' @param model Character string for model name
#' @param seeds A vector of seeds to run across (sequentially) for each
#' algorithm.
#' @param inits A list of lists of randomly drawn independent starting
#' values, of the same length as seeds
#' @param Nout The number of resulting iterations after thinning and warmup
#' @param Nthin The number to thin, defaults to 1.
#' @param lambda A vector of integrated time for HMC (eps times L)
#' @param delta A vector of target acceptance rates. Defaults to 0.8.
#' @param sink.console Whether to sink console output to trash file to cleanup
#' console output. Defaults to TRUE. Makes it easier to see the progress on
#' the console.
#' @param max_treedepth The maximum times the NUTS algorithm can double
#' before exiting that iteration. Default is 10 in Stan and this function.
#' @param metric A vector of metrics to run across for HMC and NUTS. Must
#' be one of c("unit_e", "diag_e", "dense_e").
#' @return A list of two data frames. adapt is the adaptive results from
#' Stan, and perf is the performance metrics for each run.
run.chains <- function(model, seeds, Nout, Nthin=1, lambda, delta=.8,
                       metric='diag_e', data, inits, params.jags, max_treedepth=10,
                       sink.console=TRUE){
  if(Nthin!=1) stop('this probably breaks if Nthin!=1')
  model.jags <- paste0(model, '.jags')
  model.stan <- paste0(model, '.stan')
  Niter <- 2*Nout*Nthin
  Nwarmup <- Niter/2
  ind.warmup <- 1:Nwarmup              # index of samples, excluding warmup
  ind.samples <- (Nwarmup+1):Niter     # index of warmup samples
  perf.list <- list()
  adapt.list <- list()
  k <- 1
  if(sink.console){
    sink(file='trash.txt', append=FALSE, type='output')
    on.exit(sink())
  }
  ## Precompile Stan model so it isn't done repeatedly and isn't in the
  ## timings
  model.stan <- stan(file=model.stan, data=data, iter=100, par=params.jags,
                     warmup=50, chains=1, thin=1, algorithm='NUTS',
                     init=list(inits[[1]]), seed=1, verbose=FALSE,
                     control=list(adapt_engaged=FALSE))
  for(seed in seeds){
    message(paste('==== Starting seed',seed, 'at', Sys.time()))
    inits.seed <- list(inits[[which(seed==seeds)]])
    set.seed(seed)
    ## Precompile JAGS model so not in timings
    ## message('Starting JAGS model')
    time.jags <- as.vector(system.time(fit.jags <-
    jags(data=data, parameters.to.save=params.jags, inits=inits.seed,
         model.file=model.jags, n.chains=1, DIC=FALSE,
         n.iter=Niter, n.burnin=Nwarmup, n.thin=1)))[3]
    saveRDS(fit.jags, file=paste('fits/jags', seed,'.RDS', sep='_'))
    sims.jags <- fit.jags$BUGSoutput$sims.array
    perf.jags <- data.frame(rstan::monitor(sims=sims.jags, warmup=0, print=FALSE, probs=.5))
    Rhat.jags <- with(perf.jags, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
    perf.list[[k]] <-
      data.frame(platform='jags', seed=seed, Npar=dim(sims.jags)[3], delta.target=.5,
                 time=time.jags,
                 minESS=min(perf.jags$n_eff), medianESS=median(perf.jags$n_eff),
                 Nsims=dim(sims.jags)[1],
                 minESS.coda=min(coda::effectiveSize(x=sims.jags[,1,])),
                 Rhat.jags)
    k <- k+1
    rm(sims.jags, perf.jags)
    ## Use adaptation for eps and diagonal covariances, but remove those
    ## samples and time it took
    ## message('Starting stan.nuts models')
    for(idelta in delta){
      for(imetric in metric){
        time.stan <- as.vector(system.time(fit.stan.nuts <-
          stan(fit=model.stan, data=data, iter=Niter,
               warmup=Nwarmup, chains=1, thin=Nthin, algorithm='NUTS',
               init=inits.seed, seed=seed, par=params.jags,
               control=list(adapt_engaged=TRUE, adapt_delta=idelta,
                 metric=imetric, max_treedepth=max_treedepth))))[3]
        saveRDS(fit.stan.nuts, file=paste('fits/stan_nuts', metric, idelta, seed,'.RDS', sep='_'))
        sims.stan.nuts <- extract(fit.stan.nuts, permuted=FALSE)
        perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
        Rhat.stan.nuts <- with(perf.stan.nuts, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
        adapt.nuts <- as.data.frame(get_sampler_params(fit.stan.nuts, inc_warmup=FALSE))
        adapt.list[[k]] <-
          data.frame(platform='stan.nuts', seed=seed,
                     Npar=dim(sims.stan.nuts)[3]-1,
                     Nsims=dim(sims.stan.nuts)[1],
                     delta.mean=mean(adapt.nuts$accept_stat__),
                     delta.target=idelta,
                     eps.final=tail(adapt.nuts$stepsize__,1),
                     max_treedepths=sum(adapt.nuts$treedepth__>max_treedepth),
                     ndivergent=sum(adapt.nuts$divergent__),
                     nsteps.mean=mean(adapt.nuts$n_leapfrog__),
                     nsteps.median=median(adapt.nuts$n_leapfrog__),
                     nsteps.sd=sd(adapt.nuts$n_leapfrog__),
                     metric=imetric)
        perf.list[[k]] <-
          data.frame(platform='stan.nuts',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=tail(adapt.nuts$stepsize__,1),
                     Npar=dim(sims.stan.nuts)[3]-1,
                     time=time.stan,
                     minESS=min(perf.stan.nuts$n_eff),
                     medianESS=median(perf.stan.nuts$n_eff),
                     Nsims=dim(sims.stan.nuts)[1],
                     minESS.coda=min(effectiveSize(as.data.frame(sims.stan.nuts[,1,]))),
                     Rhat.stan.nuts)
        k <- k+1
      }}
    rm(fit.stan.nuts, sims.stan.nuts, perf.stan.nuts, adapt.nuts)
### Currently turned off: feature to explore static HMC.
    ## if(!is.null(lambda)){
    ##   for(ilambda in lambda){
    ##     message(paste0('Starting stan.hmc',ilambda,' models'))
    ##     for(idelta in delta){
    ##       for(imetric in metric){
    ##         fit.stan.hmc <-
    ##           stan(fit=model.stan, data=data, iter=Niter,
    ##                warmup=Nwarmup, chains=1, thin=Nthin, algorithm='HMC', seed=seed, init=inits.seed,
    ##                par=params.jags, control=list(adapt_engaged=TRUE,
    ##                                   adapt_delta=idelta, metric=imetric, int_time=ilambda))
    ##         sims.stan.hmc <- extract(fit.stan.hmc, permuted=FALSE)
    ##         perf.stan.hmc <- data.frame(rstan::monitor(sims=sims.stan.hmc, warmup=0, print=FALSE))
    ##         Rhat.stan.hmc <- with(perf.stan.hmc, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
    ##         adapt.hmc <- as.data.frame(get_sampler_params(fit.stan.hmc))
    ##         adapt.list[[k]] <-
    ##           data.frame(platform=paste0('HMC',ilambda), seed=seed,
    ##                      Npar=dim(sims.stan.hmc)[3]-1,
    ##                      Nsims=dim(sims.stan.hmc)[1],
    ##                      delta.target=idelta,
    ##                      delta.mean=mean(adapt.hmc$accept_stat__[ind.samples]),
    ##                      eps.final=tail(adapt.hmc$stepsize__,1),
    ##                      metric=imetric, lambda=ilambda,
    ##                      L=tail(adapt.hmc$stepsize__,1)*ilambda,
    ##                      nsteps.mean=tail(adapt.hmc$stepsize__,1)*ilambda)
    ##         perf.list[[k]] <-
    ##           data.frame(platform=paste0('stan.hmc',ilambda),
    ##                      seed=seed, delta.target=idelta, metric=imetric,
    ##                      Npar=dim(sims.stan.hmc)[3]-1,
    ##                      time.sampling=get_elapsed_time(fit.stan.hmc)[2],
    ##                      time.warmup=get_elapsed_time(fit.stan.hmc)[1],
    ##                      minESS=min(perf.stan.hmc$n_eff),
    ##                      medianESS=median(perf.stan.hmc$n_eff),
    ##                      Nsims=dim(sims.stan.hmc)[1],
    ##                      minESS.coda=min(effectiveSize(as.data.frame(sims.stan.hmc[,1,]))),
    ##                      Rhat.stan.hmc)
    ##         k <- k+1
    ##         rm(fit.stan.hmc, sims.stan.hmc, perf.stan.hmc, adapt.hmc)
    ##       }
    ##     }
    ##   }
    ## }
  }
  perf <- do.call(rbind.fill, perf.list)
  perf <- within(perf, {
                   time.total <- time
                   samples.per.time <- minESS/time.total
                 })
  adapt <- do.call(rbind.fill, adapt.list[!ldply(adapt.list, is.null)])
  perf$model <- adapt$model <- model
  return(invisible(list(adapt=adapt, perf=perf)))
}

#' Verify models and then run empirical tests across delta
fit.empirical <- function(model, params.jags, model.jags, inits, data, seeds,
                          delta, lambda, model.stan, Nout,  metric,
                          Nthin=1, sink.console=TRUE, ...){
    ## Now rerun across gradient of acceptance rates and compare to JAGS
    message('Starting empirical runs')
    results.empirical <-
      run.chains(model=model, seeds=seeds, Nout=Nout, lambda=lambda,
                 metric=metric, delta=delta, data=data,
                 Nthin=Nthin, inits=inits, params.jags=params.jags,
                 sink.console=sink.console, ...)
    with(results.empirical, plot.empirical.results(perf, adapt))
    write.csv(file=results.file(paste0(m, '_adapt_empirical.csv')), results.empirical$adapt)
    write.csv(file=results.file(paste0(m, '_perf_empirical.csv')), results.empirical$perf)
}

#' Make plots comparing the performance of simulated data for a model.
#'
#' @details Makes plots with Npar (model complexity) on the x-axis and a
#' variety of things on the y-axis.
#' @param perf A data frame containing model performance data, as returned
#' by run.chains
#' @param adapt A data frame containing adaptation information from Stan,
#' as returned by run.chains.
#'
#' @return Nothing. Makes plots in local folder of performance comparisons
#' and adaptation results
plot.simulated.results <- function(perf, adapt){
    model.name <- as.character(perf$model[1])
    perf.long <- melt(perf,
                      id.vars=c('model', 'platform', 'seed', 'Npar', 'Nsims'),
                      measure.vars=c('time.total', 'minESS', 'samples.per.time'))
    perf.long <- ddply(perf.long, .(platform, Npar, variable), mutate,
                       mean.value=mean(value))
    g <- ggplot(perf.long, aes(Npar, log(value), group=platform, color=platform)) +
        geom_point()+
            geom_line(data=perf.long, aes(Npar, log(mean.value))) +
                facet_wrap('variable') + ggtitle("Performance Comparison")
    ggsave(paste0('plots/', model.name, '_perf_simulated.png'), g, width=ggwidth, height=ggheight)
    adapt$pct.divergent <- with(adapt, ndivergent/Nsims)
    adapt$pct.max.treedepths <- with(adapt, max_treedepths/Nsims)
    adapt.long <- melt(adapt,
                      id.vars=c('model', 'platform', 'seed', 'Npar', 'Nsims'),
                      measure.vars=c('delta.mean', 'eps.final',
                                     'pct.divergent', 'pct.max.treedepths'))
    adapt.long <- ddply(adapt.long, .(platform, Npar, variable), mutate,
                       mean.value=mean(value))
    g <- ggplot(adapt.long, aes(Npar, value, group=platform, color=platform)) +
        geom_point() + geom_line(data=adapt.long, aes(Npar, mean.value)) +
                facet_wrap('variable') + ggtitle("Performance Comparison")
    ggsave(paste0('plots/',model.name, '_adapt_simulated.png'), g, width=ggwidth, height=ggheight)
}

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
           measure.vars=c('time', 'minESS', 'samples.per.time'))
    perf.long$seed2 <- as.factor(with(perf.long, paste(seed, metric, sep="_")))
    g <- ggplot(subset(perf.long, platform!='jags'), aes(delta.target, log(value), group=seed2, color=metric))+
      geom_line() + geom_point() + facet_grid(variable~platform, scales='free_y') + xlim(0,1)
    ggsave(paste0('plots/',model.name, '_perf_empirical.png'), g, width=ggwidth, height=ggheight)
    adapt.long <- melt(adapt, id.vars=c('platform', 'seed', 'delta.target', 'metric'),
                       measure.vars=c('eps.final', 'delta.mean', 'nsteps.mean'))
    adapt.long$seed2 <- as.factor(with(adapt.long, paste(seed, metric, sep="_")))
    g <- ggplot(adapt.long, aes(delta.target, value, group=seed2, color=metric))+
      geom_line() + facet_grid(variable~platform, scales='free_y') + xlim(0,1)
    ggsave(paste0('plots/',model.name, '_adapt_empirical.png'), g, width=ggwidth, height=ggheight)
  }


#' Verify the models are the same (coding errors) by plotting QQ plots of
#' each parameter between Stan and JAGS.
#'
#' @param sims.stan A data frame of MCMC samples from a single chain of
#' a Stan run, using extract(fit.stan).
#' @param sims.jags A data frame of MCMC samples from a single chain of a
#' JAGS run.
#' @param perf.platforms A wide data frame with platform, variable and
#' values for Rhat and n_eff, one for each variable of a chain for each
#' platform. As created by verify.model.
#' @return ggplot object invisibly. Also makes plot in folder 'plots' in
#' current working directory.
plot.model.comparisons <- function(sims.stan, sims.jags, perf.platforms=NULL){
    ## Clean up names so they match exactly
    names(sims.stan) <- gsub('\\.', '', x=names(sims.stan))
    names(sims.jags) <- gsub('\\.', '', x=names(sims.jags))
    sims.stan$lp__ <- sims.jags$deviance <- NULL
    par.names <- names(sims.jags)
    sims.jags <- sims.jags[,par.names]
    ## Massage qqplot results into long format for ggplot
    qq <- ldply(par.names, function(i){
        temp <- as.data.frame(qqplot(sims.jags[,i], sims.stan[,i], plot.it=FALSE))
        return(cbind(par=i,temp))
    })
    ## Since can be too many parameters, break them up into pages. Stolen
    ## from
    ## http://stackoverflow.com/questions/22996911/segment-facet-wrap-into-multi-page-pdf
    noVars <- length(par.names)
    noPlots <- 25
    plotSequence <- c(seq(0, noVars-1, by = noPlots), noVars)
    ## pdf('plots/model_comparison_qqplots.pdf', onefile=TRUE,
    ## width=ggwidth,height=ggheight)
    png('plots/model_comparison_qqplots%02d.png', units='in', res=500,
        width=ggwidth,height=ggheight)
    for(ii in 2:length(plotSequence)){
        start <- plotSequence[ii-1] + 1;   end <- plotSequence[ii]
        tmp <- subset(qq, par %in% par.names[start:end])
        g <- ggplot(tmp, aes(x,y))+ geom_point(alpha=.5) +
            geom_abline(slope=1, col='red') + facet_wrap('par',
        scales='free', nrow=5) + xlab('jags')+ ylab('stan')
           ## theme(axis.text.x=element_blank(), axis.text.y=element_blank())
        g <- g+ theme(text=element_text(size=7))
        print(g)
    }
    dev.off()

    if(!is.null(perf.platforms)){
        g <- ggplot(perf.platforms, aes(platform, value)) +
            facet_wrap('variable', scales='free')
        g <- g + if(nrow(perf.platforms)>50) geom_violin() else geom_point()
        ggsave('plots/model_comparison_convergence.png', g, width=9, height=5)
    }
    return(NULL)
}
#' Run models with thinning to get independent samples which are then used
#' to verify the posteriors are the same, effectively checking for bugs
#' between models before doing performance comparisons
#'
verify.models <- function(model, params.jags, inits, data, Nout, Nthin,
                          sink.console=TRUE){
  message('Starting independent runs')
  if(sink.console){
    sink(file='trash.txt', append=FALSE, type='output')
    on.exit(sink())
  }
  Niter <- 2*Nout*Nthin
  Nwarmup <- Niter/2
  model.jags <- paste0(model, '.jags')
  model.stan <- paste0(model, '.stan')
  ## Save the algorithms used by jags
  temp <- list.samplers(jags.model(file=model.jags, data=data, inits=inits, quiet=TRUE))
  write.table(table(names(temp)), file='jags.samplers.csv', sep=',', row.names=FALSE)
  fit.jags <- jags(data=data, inits=inits, param=params.jags,
                   model.file=model.jags, n.chains=1, n.burnin=Nwarmup, n.iter=Niter,
                   n.thin=Nthin)
  sims.jags <- fit.jags$BUGSoutput$sims.array
  perf.jags <- data.frame(rstan::monitor(sims=sims.jags, warmup=0, print=FALSE, probs=.5))
  fit.stan <- stan(file=model.stan, data=data, iter=Niter, chains=1,
                   warmup=Nwarmup, thin=Nthin, init=inits, pars=params.jags)
  sims.stan <- extract(fit.stan, permuted=FALSE)
  perf.stan <- data.frame(rstan::monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5))
  perf.platforms <- rbind(cbind(platform='jags',perf.jags),
                          cbind(platform='stan',perf.stan))
  perf.platforms <- melt(perf.platforms, c('Rhat', 'n_eff'), id.vars='platform')
  plot.model.comparisons(as.data.frame(sims.stan[,1,]),
                         as.data.frame(sims.jags[,1,]), perf.platforms)
  ## Save independent samples for use in intial values later.
  sims.ind <- data.frame(sims.stan[,1,])
  sims.ind <- sims.ind[, names(sims.ind) != 'lp__']
  saveRDS(sims.ind, file='sims.ind.RDS')
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

## Growth model functions
sample.vbgf <- function(ages, Linf, k,  t0, sigma.obs){
    lengths <- Linf*(1-exp(-k*(ages-t0)))
    loglengths <- log(lengths)+ rnorm(n=length(lengths), mean=0, sd=sigma.obs)
    data.frame(ages=ages, loglengths=loglengths)
}
sample.ages <- function(n.ages, t0, Ntime) {sample((t0+1):Ntime, size=n.ages, replace=FALSE)}
sample.lengths <- function(Nfish, n.ages, logLinf.mean, logLinf.sigma,
                           logk.mean, logk.sigma, sigma.obs, t0){
    Linf.vec <- exp(logLinf.mean + rnorm(n=Nfish, 0, sd=logLinf.sigma))
    k.vec <- exp(logk.mean +rnorm(n=Nfish, mean=0, sd=logk.sigma))
    dat <- ldply(1:Nfish, function(i)
        cbind(fish=i, sample.vbgf(ages=sample.ages(n.ages, t0=t0, Ntime=Ntime),
              Linf=Linf.vec[i], k=k.vec[i], sigma.obs=sigma.obs, t0=t0)))
   return( dat)
}
ss_logistic.traj <- function(r, K, num.years, sd.catch, prop.caught,
                             years.obs, sd.obs, sd.process, F, plot=TRUE){
    catches <- trajectory <- rep(0,num.years)
    u <- rnorm(n=num.years, mean=0, sd=sd.process)
    trajectory[1] <- K + u[1]
    for( yr in 2: num.years){
        Nt <- trajectory[yr-1]
        Ct <- catches[yr-1] <- Nt*(1-exp(-ifelse(Nt<K/10, 0, F)))
        trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct) * exp(u[yr]-sd.process^2/2)
        if( trajectory[yr]<=0 ) { trajectory[yr] <- NA; break}
    }
    ## generate lognormal samples
    log.pop.obs <-
        rnorm(n=length(years.obs),
              mean=log(trajectory[years.obs]),sd=sd.obs) -sd.obs^2/2
    if(plot){
        plot(1:num.years, ylim=c(0,K*1.1),y= trajectory, type='l')
        points(1:num.years, catches, type='h')
        points(years.obs, exp(log.pop.obs), col='red')
    }
    return(list(N=num.years, catches=catches, u=u, logcpue=log.pop.obs))
    ## return(trajectory)
}
