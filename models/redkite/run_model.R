## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
data <- readRDS('data.RDS')
params.jags <- c("sjuv", "ssub", "sad", "rjuv", "rad")
inits <- list(list(sjuv=.4, ssub=.7, sad=.85, rjuv=.05, rad=.1))

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
  verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
              Nout=Nout.ind, Nthin=Nthin.ind)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i)
  list(sjuv=sims.ind$sjuv[i],
    ssub=sims.ind$ssub[i],
    sad=sims.ind$sad[i],
    rjuv=sims.ind$rjuv[i],
    rad=sims.ind$rad[i]))


## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

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

