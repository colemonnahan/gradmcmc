## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
## The true values used for empirical and simulated data
data <- readRDS('data.RDS')
inits <- list(list(a=rep(3.5, len=data$K-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9,
                   fameffphi=rep(0, len=data$nfam),
                   fameffp=rep(0, len=data$nfam),
                   yeareffphi=rep(0, len=4)))
params.jags <-
    c('a', 'a1', 'b0', 'b1', 'sigmayearphi', 'sigmaphi', 'sigmap',
      'fameffphi', 'fameffp', 'yeareffphi')

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
              Nout=Nout.ind, Nthin=Nthin.ind)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i)
    list(
        a=as.numeric(sims.ind[i, grep('a\\.', x=names(sims.ind))]),
        a1=sims.ind$a1[i],
        b0=as.numeric(sims.ind[i, grep('b0\\.', x=names(sims.ind))]),
        b1=as.numeric(sims.ind[i, grep('b1\\.', x=names(sims.ind))]),
        sigmayearphi=sims.ind$sigmayearphi[i],
        sigmaphi=sims.ind$sigmaphi[i],
        sigmap=sims.ind$sigmap[i],
        fameffphi=as.numeric(sims.ind[i, grep('fameffphi\\.', x=names(sims.ind))]),
        fameffp=as.numeric(sims.ind[i, grep('fameffp\\.', x=names(sims.ind))]),
        yeareffphi=as.numeric(sims.ind[i, grep('yeareffphi\\.', x=names(sims.ind))])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

message(paste('Finished with model:', m))
setwd('../..')

