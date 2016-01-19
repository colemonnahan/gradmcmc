## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.

## Load empirical data and inits
dat <- read.csv('tuna_data.csv')
cpue <- dat$CPUE
catches <- dat$Catches
data <- list(catches=catches, logcpue=log(cpue), N=nrow(dat))
inits <- list(list(logr=log(.8), logK=log(279), iq=5, isigma2=100, itau2=100,
             u=rep(0, len=nrow(dat))))
params.jags <- c('logr', 'logK', 'isigma2', 'iq', 'itau2', 'u')

## Get independent samples from each model to make sure they are coded the
## same
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              Nout=Nout, Nout.ind=Nout.ind, Nthin.ind=Nthin.ind,
              delta=delta.vec, lambda=lambda.vec)
summary(sims.jags)
## Now loop through model sizes and run for default parameters, using JAGS
## and NUTS only.
Npar.vec <- c(5, 10, 20)
adapt.list <- perf.list <- list()
for(i in seq_along(Npar.vec)){
    Npar <- Npar.vec[i]
    ## Reproducible data since seed set inside the function
    message(paste("======== Starting Npar=", Npar))
    set.seed(115)
    source("generate_data.R")
    temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
                   seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL)
    adapt.list[[i]] <- temp$adapt
    perf.list[[i]] <- temp$perf
    ## Save them as we go in case it fails
    perf <- do.call(rbind, perf.list)
    adapt <- do.call(rbind, adapt.list)
    write.csv(x=perf, file=results.file(paste0(m,'_perf_simulated.csv')))
    write.csv(x=adapt, file=results.file(paste0(m,'_adapt_simulated.csv')))
    rm(temp)
}
message("Making plots...")
plot.simulated.results(perf, adapt)
