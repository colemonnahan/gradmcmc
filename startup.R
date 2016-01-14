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
#' @param sink Whether to sink console output to trash file to cleanup
#' console output. Defaults to TRUE. Makes it easier to see the progress on
#' the console.
#' @return A list of two lists. adapt.list is the adaptive results from
#' Stan, and perf.list is the performance metrics for each run.
run.chains <- function(model.name, seeds, Nout, L,
                       model.jags, data.jags, inits.jags, params.jags,
                       model.stan, data.stan, inits.stan, params.stan,
                       n.burnin, n.thin, sink.console=TRUE){
    perf.list <- list()
    adapt.list <- list()
    k <- 1
    if(sink.console){
        sink(file='trash.txt', append=FALSE, type='output')
        on.exit(sink())
    }
    for(seed in seeds){
        ## Now run single long chains without thinning and timing to get
        ## performance (minESS/time) for each of the methods
        ## Run a long one to ensure good tuning
        message(paste('==== Starting seed',seed, 'at', Sys.time()))
        set.seed(seed)
        message('Starting JAGS model')
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


plot.model.results <- function(perf, adapt){
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
