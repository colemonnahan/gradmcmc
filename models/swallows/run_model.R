## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
data <- readRDS('data.RDS')
inits <- readRDS('inits.RDS')
params.jags <- names(inits[[1]])

## Get independent samples from each model to make sure they are coded the
## same
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta.vec, metric=metric,
              Nout=Nout, Nout.ind=Nout.ind, Nthin.ind=Nthin.ind,
              verify=FALSE)

## ## Now loop through model sizes and run for default parameters, using JAGS
## ## and NUTS only.
## adapt.list <- perf.list <- list()
## for(i in seq_along(Npar.vec)){
##     Npar <- Npar.vec[i]
##     ## Reproducible data since seed set inside the function
##     message(paste("======== Starting Npar=", Npar))
##     set.seed(115)
##     source("generate_data.R")
##     temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
##                    seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL)
##     adapt.list[[i]] <- temp$adapt
##     perf.list[[i]] <- temp$perf
##     ## Save them as we go in case it fails
##     perf <- do.call(rbind, perf.list)
##     adapt <- do.call(rbind, adapt.list)
##     plot.simulated.results(perf, adapt)
##     write.csv(x=perf, file=results.file(paste0(m,'_perf_simulated.csv')))
##     write.csv(x=adapt, file=results.file(paste0(m,'_adapt_simulated.csv')))
##     rm(temp)
## }
message(paste('Finished with model:', m))
setwd('../..')

