### ------------------------------------------------------------
## The purpose of this file is to prepare the models for running later. May
## need to refresh/create data or not depending on the model. The models
## are loaded into the workspace and used later on.
## Generate data
logLinf.mean <- log(50)
logk.mean <- log(.1)
t0 <- 5
logLinf.sigma <- .1
logk.sigma <- .2
Ntime <- 100
sigma.obs <- .08
sample.vbgf <- function(ages, Linf, k,  t0=5){
    lengths <- Linf*(1-exp(-k*(ages-t0)))
    loglengths <- log(lengths)+ rnorm(n=length(lengths), mean=0, sd=sigma.obs)
    data.frame(ages=ages, loglengths=loglengths)
}
sample.ages <- function(n.ages) {sample(6:Ntime, size=n.ages, replace=FALSE)}
sample.lengths <- function(Nfish, n.ages){
    Linf.vec <- exp(logLinf.mean + rnorm(n=Nfish, 0, sd=logLinf.sigma))
    k.vec <- exp(logk.mean +rnorm(n=Nfish, mean=0, sd=logk.sigma))
    dat <- ldply(1:Nfish, function(i) cbind(fish=i, sample.vbgf(ages=sample.ages(n.ages), Linf=Linf.vec[i], k=k.vec[i])))
    ## saveRDS(dat, paste0('growth_data_',Nfish,'.RDS'))
    dat
}
## set.seed(5)
## sigma.obs <- .0001
## dat.plot <- sample.lengths(Nfish=50, n.ages=5)
## g <- ggplot(dat.plot, aes(ages, lengths, group=fish)) +
##     geom_point(alpha=.5)
## ggsave('plots/simulated_growth.png', g, width=9, height=5)
## plot(ddply(dat.plot, .(ages), summarize, CV=sd(lengths)/mean(lengths)))
## sigma.obs <- .08
## set.seed(5)
## dat.plot <- sample.lengths(Nfish=50, n.ages=5)
## g <- ggplot(dat.plot, aes(ages, lengths, group=fish))+
##     geom_point(alpha=.5)
## ggsave('plots/observed_growth.png', g, width=9, height=5)
dat <- sample.lengths(Nfish=Nfish, n.ages=5)
data <- list(Nfish=Nfish, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
init <- list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf=rep(logLinf.mean, len=Nfish),
                  logk=rep(logk.mean, len=Nfish))
### ------------------------------------------------------------
message("Loading growth models into the workspace")
## JAGS models
data.jags <- data
params.jags <-
    c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma",
      "sigma_obs", "logk", "logLinf")
inits.jags <- list(init)
model.jags <- function(){
    ## hyperpriors
    logLinf_mean~dnorm(log(50),.5)
    logLinf_sigma~dunif(0,.5)
    logk_mean~dnorm(log(.1),.1)
    logk_sigma~dunif(0,.5)
    ## priors
    sigma_obs~dunif(0,.5)
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Nfish){
        logLinf[i]~dnorm(logLinf_mean, pow(logLinf_sigma, -2))
        logk[i]~dnorm(logk_mean, pow(logk_sigma, -2))
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        Linf[i] <- exp(logLinf[fish[i]])
        k[i] <- exp(logk[fish[i]])
        ypred[i] <- Linf[i]*(1-exp(-k[i]*(ages[i]-5)))
        ## ypred[i] <- logLinf[fish[i]]+
        ##     log( (1-exp(-exp(logk[fish[i]])*ages[i])))
        ## Likelihood of data
        loglengths[i]~dnorm(log(ypred[i]), pow(sigma_obs, -2))
    }
}
## temp <- jags(data=data.jags, inits=list(init), param=params.jags, model.file=model.jags,
##      n.chains=1, n.burnin=1000, n.iter=5000, n.thin=1)
## End of JAGS
### ------------------------------------------------------------

### ------------------------------------------------------------
## Stan models
data.stan <- data
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
inits.stan <- list(init)
model.stan <- stan(file='growth.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)
## temp <- stan(fit=model.stan, data=data.stan, iter=5000, chains=1,
##                    warmup=1000, thin=1, init=inits.stan)
rm(dat, data, init,  logk.mean, logk.sigma, logLinf.mean, logLinf.sigma,
   sample.ages, sample.lengths, sample.vbgf, sigma.obs, t0)
## End of Stan
### ------------------------------------------------------------
message("Finished loading growth models")


### ------------------------------------------------------------
## OLD CODE for testing
## model.R <- function(Linf_mean, Linf_sigma, k_mean, k_sigma, sigma_obs,
##                     Linf, k){
##     ages <- data$ages; fish <- data$fish; Nobs=data$Nobs
##     loglengths <- data$loglengths
##     NLL <- 0
##     ## hyperpriors
##     NLL <- NLL-dunif(x=Linf_mean, 30,70, log=TRUE)
##     NLL <- NLL-dunif(x=Linf_sigma, 0,0.05, log=TRUE)
##     NLL <- NLL-dunif(x=k_mean, 0, .5, log=TRUE)
##     NLL <- NLL-dunif(x=k_sigma, 0,0.05, log=TRUE)
##     ## priors
##     NLL <- NLL-dunif(sigma_obs, 0,1)
##     ## Loop through the hyperparameters (on group) and calculate
##     ## probabilities
##     NLL <- NLL- sum(dnorm(Linf, Linf_mean, Linf_sigma, log=TRUE))
##     NLL <- NLL- sum(dnorm(k, k_mean, k_sigma, log=TRUE))
##     ypred <- rep(NA, len=Nobs)
##     for(i in 1:Nobs){
##         ypred[i] <- log(Linf[fish[i]])+log( (1-exp(-k[fish[i]]*ages[i])) )
##         ## Likelihood of data
##         NLL <- NLL-dnorm(loglengths[i], 1.5, sigma_obs, log=TRUE)
##     }
##     NLL
## }
## model.R(Linf.mean, Linf.sigma, k.mean, k.sigma, sigma.obs, Linf, k)
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
##   dyn.load(TMB::dynlib('growth'))
##   model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb.est, DLL='growth')
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
## growth.nuts.perf.list <- sfLapply(1:nrow(params), function(i){
##   dyn.load(TMB::dynlib("growth"))
##   model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='growth')
##    covar.temp <- if(params$covar[i]) covar.tmb else NULL
##    run_nuts(obj=model.tmb, nsim=nsim, seed=params$seed[i],
##             covar=covar.temp, Madapt=Madapt, delta=params$delta[i])
## })
## growth.nuts.perf <- do.call(rbind, growth.nuts.perf.list)
## saveRDS(growth.nuts.perf, results.file('growth.nuts.perf.RDS'))
## growth.nuts.perf$acceptance <- NULL
## growth.nuts.perf.long <-
##     reshape2::melt(growth.nuts.perf, c('covar', 'tuning', 'algorithm', 'seed'))
## ggplot(growth.nuts.perf.long, aes(tuning, value, group=covar, color=covar)) +
##     geom_line() + facet_wrap('variable', nrow=3, scales='free')
## ggsave(plots.file('growth_nuts_perf.png'), width=8, height=6)
## ### Old way of using plyr, but can't get parallel to work
## ## growth.nuts.perf <-
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
## message("Starting TMB models for growth")
## ## Try combinations of TMB
## tmb.nuts.eps <- 0.05
## tmb.nuts.covar.eps <- .7
## tmb.hmc.eps <- 0.05
## tmb.hmc.covar.eps <- .7
## model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='growth')
## time.tmb.nuts <-
##     as.vector(system.time(
##         tmb.nuts <-
##             mcmc(obj=model.tmb, nsim=n.iter1, eps=tmb.nuts.eps,
##                  max_doubling=10, Madapt=n.burnin1, delta=.65,
##                  algorithm='NUTS', diag=FALSE)))[3]
## ## tmb.nuts <- tmb.nuts[-(1:n.burnin1),]
## tmb.nuts$LP <- -apply(tmb.nuts, 1, model.tmb$fn)
## perf.list[['tmb.nuts']] <-
##     data.frame(model='growth', platform='tmb.nuts',time=time.tmb.nuts,
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
##     data.frame(model='growth', platform='tmb.nuts.covar',time=time.tmb.nuts.covar,
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
##     data.frame(model='growth', platform='tmb.hmc1',time=time.tmb.hmc1,
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
##     data.frame(model='growth', platform='tmb.hmc1.covar',time=time.tmb.hmc1.covar,
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
##     data.frame(model='growth', platform='tmb.hmc10',time=time.tmb.hmc10,
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
##     data.frame(model='growth', platform='tmb.hmc10.covar',time=time.tmb.hmc10.covar,
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
##     data.frame(model='growth', platform='tmb.rwm.covar',time=time.tmb.rwm.covar,
##                minESS=min(effectiveSize(tmb.rwm.covar)),
##                medianESS=median(effectiveSize(tmb.rwm.covar)),
##                n=nrow(tmb.rwm.covar))
## yy <- data.frame(extract(temp, permuted=FALSE)[,1,])
## par(mfrow=c(2,3))
## myre <- function(x,y) (x-y)/x
## myxlim <- c(-.5,.5)
## hist(myre(yy$logLinf_mean,logLinf.mean), xlim=myxlim)
## hist(myre(yy$logLinf_sigma,logLinf.sigma), xlim=myxlim)
## hist(myre(yy$logk_mean,logk.mean), xlim=myxlim)
## hist(myre(yy$logk_sigma,logk.sigma), xlim=myxlim)
## hist(myre(yy$sigma_obs,sigma.obs), xlim=myxlim)
## plot(data$ages, data$loglengths, pch='.')
## a <- 5:Ntime
## lines(a, log(exp(logLinf.mean)*(1-exp(-exp(logk.mean)*(a-5)))))
## lines(a, log(exp(mean(yy$logLinf_mean))*(1-exp(-exp(mean(yy$logk_mean))*(a-5)))))
## acf(xx)
## acf(yy[,1:6])
## par(mfrow=c(2,3))
## qqplot(yy$logLinf_mean, xx$logLinf_mean); abline(0,1)
## qqplot(yy$logLinf_sigma, xx$logLinf_sigma); abline(0,1)
## qqplot(yy$logk_mean, xx$logk_mean); abline(0,1)
## qqplot(yy$logk_sigma, xx$logk_sigma); abline(0,1)
## qqplot(yy$sigma_obs, xx$sigma_obs); abline(0,1)
## xx <- data.frame(temp$BUGSoutput$sims.list)
## ## acf(xx)
## par(mfrow=c(2,3))
## myre <- function(x,y) (x-y)/x
## myxlim <- c(-.5,.5)
## hist(myre(xx$logLinf_mean,logLinf.mean), xlim=myxlim)
## hist(myre(xx$logLinf_sigma,logLinf.sigma), xlim=myxlim)
## hist(myre(xx$logk_mean,logk.mean), xlim=myxlim)
## hist(myre(xx$logk_sigma,logk.sigma), xlim=myxlim)
## hist(myre(xx$sigma_obs,sigma.obs), xlim=myxlim)
## plot(data$ages, data$loglengths, pch='.')
## a <- 5:Ntime
## lines(a, log(exp(logLinf.mean)*(1-exp(-exp(logk.mean)*(a-5)))))
## lines(a, log(exp(mean(xx$logLinf_mean))*(1-exp(-exp(mean(xx$logk_mean))*(a-5)))))

