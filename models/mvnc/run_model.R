## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
Npar <- 2
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
    verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
              Nout=Nout.ind, Nthin=Nthin.ind)

## Use independent draws from the verify.model output to use for initial
## values. Each model has a different way to format the inits data.
sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i) list(mu=as.numeric(sims.ind[i,])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

## Now loop through model sizes and run for default parameters, using JAGS
## and NUTS only.
adapt.list <- perf.list <- list()
k <- 1
for(i in seq_along(Npar.vec)){
  Npar <- Npar.vec[i]
  message(paste("============ Starting Npar=", Npar))
  for(j in seq_along(cor.vec)){
    cor.temp <- cor.vec[j]
    message(paste("============ Starting cor=", cor.temp))
    ## Reproducible data since seed set inside the function
    set.seed(115)
    source("generate_data.R")
    temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
                   seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL, delta=delta)
    adapt.list[[k]] <- cbind(temp$adapt, cor=cor.temp)
    perf.list[[k]] <- cbind(temp$perf, cor=cor.temp)
    ## Save them as we go in case it fails
    perf <- do.call(rbind, perf.list)
    adapt <- do.call(rbind, adapt.list)
    plot.simulated.results(perf, adapt)
    write.csv(x=perf, file=results.file(paste0(m,'_perf_simulated.csv')))
    write.csv(x=adapt, file=results.file(paste0(m,'_adapt_simulated.csv')))
    rm(temp)
    k <- k+1
}
}
message(paste('Finished with model:', m))

setwd('../..')
