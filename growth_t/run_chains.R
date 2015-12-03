## ehook: An analysis of the Hamley and Skud (1978) data on the effect of
## hook

setwd('growth_t')
Nout <- 5000
stan.burnin <- 2000
jags.burnin <- 2000
n.thin <- 5
Nfish <- 100                           # this sets the dimension of the model
perf.list <- list()
adapt.list <- list()

for(Nfish in c(10,100, 500, 1000, 2000)){
### ------------------------------------------------------------
source("load_models.R")
## Now run single long chains without thinning and timing to get
## performance (minESS/time) for each of the methods
set.seed(1355)
## Run a long one to ensure good tuning
temp <- jags(data=data.jags, parameters.to.save=params.jags, inits=inits.jags,
         model.file=model.jags, n.chains=1, DIC=FALSE,
         n.iter=jags.burnin+10, n.burnin=jags.burnin, n.thin=1)
## Rerun with those tunings to get efficiency
time.jags <-
    as.vector(system.time(results.jags <- update(temp, n.iter=n.thin*Nout, n.thin=n.thin))[3])
sims.jags <- results.jags$BUGSoutput$sims.array
## saveRDS(sims.jags, file=paste0('results/sims.jags.', Nfish, '.RDS'))
perf.jags <- data.frame(rstan::monitor(sims=sims.jags, warmup=0, print=FALSE))
perf.list[[paste0('jags.',Nfish)]] <-
    data.frame(model='growth_t', platform='jags',time=time.jags,
               minESS=min(perf.jags$n_eff), Nfish=Nfish,
               medianESS=median(perf.jags$n_eff),
               n=dim(sims.jags)[1])
rm(temp, sims.jags, perf.jags)
## Use adaptation for eps and diagonal covariances, but remove those
## samples and time it took
message('Starting stan.nuts model')
results.stan.nuts <-
    stan(fit=model.stan, data=data.stan, iter=n.thin*Nout+stan.burnin,
         warmup=stan.burnin, chains=1, thin=n.thin, algorithm='NUTS',
         init=inits.stan,
         control=list(adapt_engaged=TRUE))
adapt.list[[paste0('stan.nuts.',Nfish)]] <- as.data.frame(get_sampler_params(results.stan.nuts))
time.stan.nuts <- get_elapsed_time(results.stan.nuts)[2]
sims.stan.nuts <- extract(results.stan.nuts, permuted=FALSE)
## saveRDS(sims.stan.nuts, file=paste0('results/sims.stan.nuts.', Nfish, '.RDS'))
perf.stan.nuts <- data.frame(rstan::monitor(sims=sims.stan.nuts, warmup=0, print=FALSE))
perf.list[[paste0('stan.nuts.',Nfish)]] <-
    data.frame(model='growth_t', platform='stan.nuts',time=time.stan.nuts,
               minESS=min(perf.stan.nuts$n_eff), Nfish=Nfish,
               medianESS=median(perf.stan.nuts$n_eff),
               n=dim(sims.stan.nuts)[1])
rm(results.stan.nuts, sims.stan.nuts, perf.stan.nuts)
message('Starting stan.hmc1 model')
results.stan.hmc1 <-
    stan(fit=model.stan, data=data.stan, iter=n.thin*Nout+stan.burnin,
         warmup=stan.burnin, chains=1, thin=n.thin, algorithm='HMC',
         init=inits.stan, control=list(adapt_engaged=TRUE, int_time=1))
adapt.list[[paste0('stan.hmc1.',Nfish)]] <- as.data.frame(get_sampler_params(results.stan.hmc1))
time.stan.hmc1 <- get_elapsed_time(results.stan.hmc1)[2]
sims.stan.hmc1 <- extract(results.stan.hmc1, permuted=FALSE)
## saveRDS(sims.stan.hmc1, file=paste0('results/sims.stan.hmc1.', Nfish, '.RDS'))
perf.stan.hmc1 <- data.frame(rstan::monitor(sims=sims.stan.hmc1, warmup=0, print=FALSE))
perf.list[[paste0('stan.hmc1.',Nfish)]] <-
    data.frame(model='growth_t', platform='stan.hmc1',time=time.stan.hmc1,
               minESS=min(perf.stan.hmc1$n_eff), Nfish=Nfish,
               medianESS=median(perf.stan.hmc1$n_eff),
               n=dim(sims.stan.hmc1)[1])
rm(results.stan.hmc1, sims.stan.hmc1, perf.stan.hmc1)
message('Starting stan.hmc10 model')
results.stan.hmc10 <-
    stan(fit=model.stan, data=data.stan, iter=n.thin*Nout+stan.burnin,
         warmup=stan.burnin, chains=1, thin=n.thin, algorithm='HMC',
         init=inits.stan, control=list(adapt_engaged=TRUE, int_time=10))
adapt.list[[paste0('stan.hmc10.',Nfish)]] <- as.data.frame(get_sampler_params(results.stan.hmc10))
time.stan.hmc10 <- get_elapsed_time(results.stan.hmc10)[2]
sims.stan.hmc10 <- extract(results.stan.hmc10, permuted=FALSE)
## saveRDS(sims.stan.hmc10, file=paste0('results/sims.stan.hmc10.', Nfish, '.RDS'))
perf.stan.hmc10 <- data.frame(rstan::monitor(sims=sims.stan.hmc10, warmup=0, print=FALSE))
perf.list[[paste0('stan.hmc10.',Nfish)]] <-
    data.frame(model='growth_t', platform='stan.hmc10',time=time.stan.hmc10,
               minESS=min(perf.stan.hmc10$n_eff), Nfish=Nfish,
               medianESS=median(perf.stan.hmc10$n_eff),
               n=dim(sims.stan.hmc10)[1])
rm(results.stan.hmc10, sims.stan.hmc10, perf.stan.hmc10)
saveRDS(perf.list, 'results/perf.list.RDS')
saveRDS(adapt.list, 'results/adapt.list.RDS')
}

message("Making  plots")
perf.list <- readRDS('results/perf.list.RDS')
perf <- do.call(rbind, perf.list)
perf$samples.per.time <- perf$minESS/perf$time
perf$pct.ess <- 100*perf$minESS/perf$n
perf$N.parameters <- with(perf, 5+2*Nfish)
perf[,c('platform', 'n')]
ggplot(perf, aes(N.parameters, time, group=platform, color=platform))+geom_line()+geom_point()
ggsave('plots/growth_t_perf_time.png', width=12, height=5)
ggplot(perf, aes(N.parameters, pct.ess, group=platform, color=platform))+geom_line()+geom_point()
ggsave('plots/growth_t_perf_pctess.png', width=12, height=5)
ggplot(perf, aes(N.parameters, medianESS, group=platform, color=platform))+geom_line()+geom_point()
ggsave('plots/growth_t_perf_medianESS.png', width=12, height=5)
ggplot(perf, aes(N.parameters, log(samples.per.time), group=platform, color=platform))+geom_line()+geom_point()
ggsave('plots/growth_t_perf_perf.png', width=12, height=5)


## ## Make sure Stan is adapting well enough
## message("Making cumulative minESS plots")
## x.jags <- cbind('software'='jags', cum.minESS(sims.jags))
## x.stan.nuts <- cbind('software'='stan.nuts', cum.minESS(sims.stan.nuts))
## x.stan.hmc1 <- cbind('software'='stan.hmc1', cum.minESS(sims.stan.hmc1))
## x.stan.hmc10 <- cbind('software'='stan.hmc10', cum.minESS(sims.stan.hmc10))
## x <- rbind(x.jags, x.stan.nuts, x.stan.hmc1, x.stan.hmc10)
## ggplot(x, aes(iteration, pct.ess, group=software, color=software))+geom_line()
## ggsave(paste0('plots/growth_t_cum_miness_',Nfish,'.png'), width=7,
## height=5)
## par(mfrow=c(3,4), mar=c(3,3,.5,.5), oma=c(0,0,0,0))
## myylim <- c(0,1)
## plot(adapt.list[['stan.nuts.10']]$stepsize__, ylim=myylim)
## mtext
## plot(adapt.list[['stan.nuts.100']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.nuts.500']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.nuts.1000']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.nuts.10']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.nuts.100']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.nuts.500']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.nuts.1000']]$accept_stat__, ylim=myylim)
## myylim <- c(0,100)
## plot(adapt.list[['stan.nuts.10']]$n_leapfrog__, ylim=myylim)
## plot(adapt.list[['stan.nuts.100']]$n_leapfrog__, ylim=myylim)
## plot(adapt.list[['stan.nuts.500']]$n_leapfrog__, ylim=myylim)
## plot(adapt.list[['stan.nuts.1000']]$n_leapfrog__, ylim=myylim)

## par(mfrow=c(2,4), mar=c(3,3,.5,.5), oma=c(0,0,0,0))
## myylim <- c(0,1)
## plot(adapt.list[['stan.hmc1.10']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.100']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.500']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.1000']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.10']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.100']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.500']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc1.1000']]$accept_stat__, ylim=myylim)

## par(mfrow=c(2,4), mar=c(3,3,.5,.5), oma=c(0,0,0,0))
## myylim <- c(0,1)
## plot(adapt.list[['stan.hmc10.10']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.100']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.500']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.1000']]$stepsize__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.10']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.100']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.500']]$accept_stat__, ylim=myylim)
## plot(adapt.list[['stan.hmc10.1000']]$accept_stat__, ylim=myylim)


## test <- do.call(rbind, adapt.list[[-grep('nuts',names(adapt.list))]])

## ### ------------------------------------------------------------
## ## Run a long chain with thinning to get independent samples to make sure
## ## the models are matching.
## ## JAGS
## stop('do not run ind samples')
## n.out.ind <- 2000
## n.thin.ind <- 1000
## n.chains.ind <- 4
## n.iter.ind <- 1.25*n.out.ind*n.thin.ind
## n.burnin.ind <- min(2000,floor(.1*n.iter.ind))
## results.jags.ind <-
##     jags(data=data.jags, parameters.to.save=params.jags,
##          model.file=model.jags, n.chains=n.chains.ind,
##          n.iter=n.iter.ind, n.burnin=n.burnin.ind, n.thin=n.thin.ind)
## saveRDS(results.jags.ind, file='results/results.jags.ind.RDS')
## ## Stan
## results.stan.ind <-
##     stan(fit=model.stan, data=data.stan, iter=n.iter.ind, warmup=n.burnin.ind,
##          chains=n.chains.ind, thin=n.thin.ind, algorithm='NUTS',
##          init=rep(list(inits), n.chains.ind),
##          control=list(adapt_engaged=TRUE))
## saveRDS(results.stan.ind, file='results/results.stan.ind.RDS')
## ## TMB
## sfStop()
## sfInit( parallel=TRUE, cpus=n.chains.ind )
## sfExport("data.tmb", "inits.tmb.est", "n.thin.ind", "n.iter.ind",
##          "n.burnin.ind")
## temp <- sfLapply(1:n.chains.ind, function(i){
##   dyn.load(TMB::dynlib('growth_t'))
##   model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb.est, DLL='growth_t')
##   set.seed(i)
##   x <- TMB::mcmc(obj=model.tmb, nsim=n.iter.ind, eps=NULL, max_doubling=6,
##             Madapt=n.burnin.ind, delta=.65, algorithm='NUTS', diag=FALSE)
##   ## par(mfrow=c(1,3));with(x, {plot(Hbar); plot(epsbar); plot(epsvec)})
##   ## plot(x$par$beta2)
##   ## discard warmup and then thin
##   x <- x[-(1:n.burnin.ind),]
##   x <- x[seq(1, nrow(x), by=n.thin.ind),]
##   x$LP <- -apply(x, 1, model.tmb$fn)
##   cbind(chain=i, x)
## })
## results.tmb.ind <- do.call(rbind, temp)
## saveRDS(results.tmb.ind, file='results/results.tmb.ind.RDS')
## ## End of independent sampling
## ### ------------------------------------------------------------


### quick check that the adaption is going well and long enough
## test <- llply(delta.vec, function(x)
##     run_nuts(obj=model.tmb, nsim=nsim, seed=seed, covar=cov,
##              Madapt=Madapt, delta=x, diag=TRUE, inits=inits.tmb))
## test2 <- data.frame(do.call(cbind, llply(test, function(x) x$epsbar)))
## names(test2) <- delta.vec
## test2$iteration <- 1:nrow(test2)
## ggplot(reshape2::melt(test2, c('iteration')),
##        aes(iteration, value, group=variable, color=variable)) +
##     geom_line()




## temp2 <- lapply(temp, function(x) x[-(1:warmup),])
## temp3 <- lapply(temp2, function(x) x[seq(1, nrow(x), by=thin),])

## nlls <- do.call(rbind, lapply(1:n.chains.ind, function(x) data.frame(iteration=1:nrow(temp3[[x]]),seed=x, nll=temp3[[x]]$nll )))
## ggplot(nlls, aes(iteration, nll, group=seed, color=factor(seed)))+geom_line()
## tmb.ind.ess <- data.frame(do.call(rbind, lapply(temp3, effectiveSize)))
## tmb.ind.ess$seed <- 1:n.chains.ind
## tmb.ind.ess.long <- reshape2::melt(tmb.ind.ess, 'seed', value.name='ess')
## ggplot(tmb.ind.ess.long, aes(variable, ess, color=factor(seed)))+geom_point()

## results.tmb.ind <- do.call(rbind, temp3)
## saveRDS(results.tmb.ind, file='results/results.tmb.ind.RDS')

## delta.vec <- seq(.1,.9, len=n.chains.ind*2)
## seeds.vec <- 1
## nsim <- 10000
## Madapt <- min(200, .1*nsim)
## params <- data.frame(expand.grid(delta=delta.vec, covar=c(TRUE, FALSE), seed=seeds.vec))
## sfExport("params", "Madapt", "nsim", "covar.tmb", "run_nuts", "data.tmb",
##          "inits.tmb", "run_hmc")
## growth_t.nuts.perf.list <- sfLapply(1:nrow(params), function(i){
##   dyn.load(TMB::dynlib("growth_t"))
##   model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='growth_t')
##    covar.temp <- if(params$covar[i]) covar.tmb else NULL
##    run_nuts(obj=model.tmb, nsim=nsim, seed=params$seed[i],
##             covar=covar.temp, Madapt=Madapt, delta=params$delta[i])
## })
## growth_t.nuts.perf <- do.call(rbind, growth_t.nuts.perf.list)
## saveRDS(growth_t.nuts.perf, results.file('growth_t.nuts.perf.RDS'))
## growth_t.nuts.perf$acceptance <- NULL
## growth_t.nuts.perf.long <-
##     reshape2::melt(growth_t.nuts.perf, c('covar', 'tuning', 'algorithm', 'seed'))
## ggplot(growth_t.nuts.perf.long, aes(tuning, value, group=covar, color=covar)) +
##     geom_line() + facet_wrap('variable', nrow=3, scales='free')
## ggsave(plots.file('growth_t_nuts_perf.png'), width=8, height=6)



## ### Old way of using plyr, but can't get parallel to work
## ## growth_t.nuts.perf <-
## ##     ldply(seeds.vec, function(seed)
## ##           ldply(covar.list, function(cov)
## ##               llply(delta.vec, function(x)
## ##                   run_nuts(obj=model.tmb, nsim=nsim, seed=seed, covar=cov,
## ##                            Madapt=Madapt, delta=x)
## ##                     )))

## ## Looking at how the DA works with NaN values for alpha
## epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
## eps <- epsvec[1] <- epsbar[1] <- .0078125
## mu <- log(10*eps)
## Hbar[1] <- .001; gamma <- 0.05; t0 <- 10; kappa <- 0.75
## for(m in 1:Madapt){
##     logalpha <- NaN
##     Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
##             (delta-min(1,exp(logalpha)))/(m+t0)
##     ## If logalpha not defined, skip this updating step and use
##     ## the last one.
##     if(is.nan(Hbar[m+1])) Hbar[m+1] <- 2*abs(Hbar[m])
##     logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
##     epsvec[m+1] <- exp(logeps)
##     logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
##     epsbar[m+1] <- exp(logepsbar)
##     eps <- epsvec[m+1]
## }
## plot(Hbar); plot(epsbar); plot(epsvec)
       ##
       ## TMB:::.find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=tmb.nuts.covar.eps)
## u <- TMB:::.sample.u(theta=theta.cur, r=r.cur, fn=fn2)
## TMB:::.buildtree(theta.cur, r=r.cur, v=1, j=0, eps=tmb.nuts.covar.eps,
##                  theta0=theta.cur, r0=r.cur, fn=fn2, gr=gr2)
## exp(TMB:::.calculate.H(theta=theta.cur,r=r.cur, fn=fn2))
## TMB:::.calculate.H(theta=theta, r=r, fn=fn)
## fn(theta.cur)-(1/2)*sum(r.cur^2)
## par <- c(-42971704273706.6, -14485844811983.9, 13861003117130.6, 14850057275244.7,
##   13870817587701.3, 13876925850253.6, 17076586228449.7, 10695806413971.9,
##   12320496167034.9, 19455141565795, 19931634722376.5, 27979336044559.6,
##   20757597588157.3, 6877188659907.08, 10778035982904.9, 10195422909598.6,
##   6089208092074.59, -1422561742193.6, 3357999633890.51, 2987375473482.92,
##   2714047127953.54, -184175140022.005, 30790845764431.3, -7572062312502.58,
##   12157785618377.1, 1208126576837.3, 4770810387536.19, 96168315204058.4,
##   758577078225.313, -4758124087735.35, -6917241097559.55,
##   -6704352992288.63, 3741578389014.47, -77368078987773.8)
## model.tmb$fn(par)
## model.tmb$report()

## message("Starting TMB models for growth_t")
## ## Try combinations of TMB
## tmb.nuts.eps <- 0.05
## tmb.nuts.covar.eps <- .7
## tmb.hmc.eps <- 0.05
## tmb.hmc.covar.eps <- .7
## model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='growth_t')
## time.tmb.nuts <-
##     as.vector(system.time(
##         tmb.nuts <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.nuts.eps,
##                  max_doubling=10, Madapt=n.burnin1, delta=.65,
##                  algorithm='NUTS', diag=FALSE)))[3]
## ## tmb.nuts <- tmb.nuts[-(1:n.burnin1),]
## tmb.nuts$LP <- -apply(tmb.nuts, 1, model.tmb$fn)
## perf.list[['tmb.nuts']] <-
##     data.frame(model='growth_t', platform='tmb.nuts',time=time.tmb.nuts,
##                minESS=min(effectiveSize(tmb.nuts)),
##                medianESS=median(effectiveSize(tmb.nuts)),
##                n=nrow(tmb.nuts))
## time.tmb.nuts.covar <-
##     as.vector(system.time(
##         tmb.nuts.covar <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.nuts.covar.eps,
##                  max_doubling=6, covar=covar.tmb, Madapt=n.burnin1,
##                  delta=.65, algorithm='NUTS', diag=FALSE)))[3]
## ## tmb.nuts.covar <- tmb.nuts.covar[-(1:n.burnin1),]
## tmb.nuts.covar$LP <- -apply(tmb.nuts.covar, 1, model.tmb$fn)
## perf.list[['tmb.nuts.covar']] <-
##     data.frame(model='growth_t', platform='tmb.nuts.covar',time=time.tmb.nuts.covar,
##                minESS=min(effectiveSize(tmb.nuts.covar)),
##                medianESS=median(effectiveSize(tmb.nuts.covar)),
##                n=nrow(tmb.nuts.covar))
## time.tmb.hmc1 <-
##     as.vector(system.time(
##         tmb.hmc1 <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.hmc.eps,
##                  Madapt=n.burnin1, delta=.65,
##                  algorithm='HMC', L=1, diag=FALSE)))[3]
## ## tmb.hmc1 <- tmb.hmc1[-(1:n.burnin1),]
## tmb.hmc1$LP <- -apply(tmb.hmc1, 1, model.tmb$fn)
## perf.list[['tmb.hmc1']] <-
##     data.frame(model='growth_t', platform='tmb.hmc1',time=time.tmb.hmc1,
##                minESS=min(effectiveSize(tmb.hmc1)),
##                medianESS=median(effectiveSize(tmb.hmc1)),
##                n=nrow(tmb.hmc1))
## time.tmb.hmc1.covar <-
##     as.vector(system.time(
##         tmb.hmc1.covar <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.hmc.covar.eps,
##                  covar=covar.tmb, Madapt=n.burnin1,
##                  delta=.65, algorithm='HMC', L=1, diag=FALSE)))[3]
## ## tmb.hmc1.covar <- tmb.hmc1.covar[-(1:n.burnin1),]
## tmb.hmc1.covar$LP <- -apply(tmb.hmc1.covar, 1, model.tmb$fn)
## perf.list[['tmb.hmc1.covar']] <-
##     data.frame(model='growth_t', platform='tmb.hmc1.covar',time=time.tmb.hmc1.covar,
##                minESS=min(effectiveSize(tmb.hmc1.covar)),
##                medianESS=median(effectiveSize(tmb.hmc1.covar)),
##                n=nrow(tmb.hmc1.covar))
## time.tmb.hmc10 <-
##     as.vector(system.time(
##         tmb.hmc10 <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.hmc.eps,
##                  Madapt=n.burnin1, delta=.65,
##                  algorithm='HMC', L=1, diag=FALSE)))[3]
## ## tmb.hmc10 <- tmb.hmc10[-(1:n.burnin1),]
## tmb.hmc10$LP <- -apply(tmb.hmc10, 1, model.tmb$fn)
## perf.list[['tmb.hmc10']] <-
##     data.frame(model='growth_t', platform='tmb.hmc10',time=time.tmb.hmc10,
##                minESS=min(effectiveSize(tmb.hmc10)),
##                medianESS=median(effectiveSize(tmb.hmc10)),
##                n=nrow(tmb.hmc10))
## time.tmb.hmc10.covar <-
##     as.vector(system.time(
##         tmb.hmc10.covar <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.hmc.covar.eps,
##                  covar=covar.tmb, Madapt=n.burnin1,
##                  delta=.65, algorithm='HMC', L=1, diag=FALSE)))[3]
## ## tmb.hmc10.covar <- tmb.hmc10.covar[-(1:n.burnin1),]
## tmb.hmc10.covar$LP <- -apply(tmb.hmc10.covar, 1, model.tmb$fn)
## perf.tmb.hmc10.covar <- time.tmb.hmc10.covar/(100*min(effectiveSize(tmb.hmc10.covar))/nrow(tmb.hmc10.covar))
## perf.list[['tmb.hmc10.covar']] <-
##     data.frame(model='growth_t', platform='tmb.hmc10.covar',time=time.tmb.hmc10.covar,
##                minESS=min(effectiveSize(tmb.hmc10.covar)),
##                medianESS=median(effectiveSize(tmb.hmc10.covar)),
##                n=nrow(tmb.hmc10.covar))
## time.tmb.rwm.covar <-
##     as.vector(system.time(
##         tmb.rwm.covar <-
##             mcmc(obj=model.tmb, nsim=n.iter1, alpha=.3,
##                  covar=covar.tmb, algorithm='RWM', diag=FALSE)))[3]
## ## tmb.rwm.covar <- tmb.rwm.covar[-(1:n.burnin1),]
## tmb.rwm.covar$LP <- -apply(tmb.rwm.covar, 1, model.tmb$fn)
## perf.tmb.rwm.covar <- time.tmb.rwm.covar/(100*min(effectiveSize(tmb.rwm.covar))/nrow(tmb.rwm.covar))
## perf.list[['tmb.rwm.covar']] <-
##     data.frame(model='growth_t', platform='tmb.rwm.covar',time=time.tmb.rwm.covar,
##                minESS=min(effectiveSize(tmb.rwm.covar)),
##                medianESS=median(effectiveSize(tmb.rwm.covar)),
##                n=nrow(tmb.rwm.covar))
