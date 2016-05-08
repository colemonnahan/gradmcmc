## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
## The true values used for empirical and simulated data
logLinf.mean <- log(50)
logk.mean <- log(.1)
t0 <- 6
logLinf.sigma <- .1
logk.sigma <- .2
Ntime <- 40
sigma.obs <- .1
Nfish <- 30
set.seed(115)
dat <- sample.lengths(Nfish=Nfish, n.ages=5, logLinf.mean=logLinf.mean,
                           logLinf.sigma=logLinf.sigma, logk.mean=logk.mean,
                           logk.sigma=logk.sigma, sigma.obs=sigma.obs, t0=t0)
g <- ggplot(dat, aes(ages, loglengths, group=fish, color=fish)) +
    geom_point(alpha=.5) + geom_line()
ggsave('plots/simulated_growth.png', g, width=9, height=5)
data <- list(Nfish=Nfish, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
inits <- list(list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf_raw=rep(0, len=Nfish),
                  logk_raw=rep(0, len=Nfish), delta=1))
params.jags <-
    c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma",
      "sigma_obs", "logk_raw","logLinf_raw", "delta")

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
  verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
                Nout=Nout.ind, Nthin=Nthin.ind, sink=sink)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i)
    list(
        logLinf_mean=sims.ind$logLinf_mean[i],
        logLinf_sigma=sims.ind$logLinf_sigma[i],
        logk_mean=sims.ind$logk_mean[i],
        logk_sigma=sims.ind$logk_sigma[i],
        sigma_obs=sims.ind$sigma_obs[i],
        delta=sims.ind$delta[i],
        logk_raw=as.numeric(sims.ind[i, grep('logk_raw', x=names(sims.ind))]),
        logLinf_raw=as.numeric(sims.ind[i, grep('logLinf_raw', x=names(sims.ind))])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

## Now loop through model sizes and run for default parameters, using JAGS
## and NUTS only.
adapt.list <- perf.list <- list()
for(i in seq_along(Npar.vec)){
    Npar <- Npar.vec[i]
    ## Reproducible data since seed set inside the function
    message(paste("======== Starting Npar=", Npar))
    set.seed(115)
    source("generate_data.R")
    temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
                   seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL, delta=delta)
    adapt.list[[i]] <- temp$adapt
    perf.list[[i]] <- temp$perf
    ## Save them as we go in case it fails
    perf <- do.call(rbind, perf.list)
    adapt <- do.call(rbind, adapt.list)
    plot.simulated.results(perf, adapt)
    write.csv(x=perf, file=results.file(paste0(m,'_perf_simulated.csv')))
    write.csv(x=adapt, file=results.file(paste0(m,'_adapt_simulated.csv')))
    rm(temp)
}
message(paste('Finished with model:', m))

setwd('../..')

