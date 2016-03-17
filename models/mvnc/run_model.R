## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
Npar <- 5
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'

## Get independent samples from each model to make sure they are coded the
## same
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta.vec, metric=metric,
              Nout=Nout, Nout.ind=Nout.ind, Nthin.ind=Nthin.ind,
              verify=TRUE)

## Now loop through model sizes and run for default parameters, using JAGS
## and NUTS only.
adapt.list <- perf.list <- list()
Npar <- 5
for(i in seq_along(cor.vec)){
    cor.temp <- cor.vec[i]
    ## Reproducible data since seed set inside the function
    message(paste("======== Starting cor=", cor.temp))
    set.seed(115)
    source("generate_data.R")
    temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
                   seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL)
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
