## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
## The true values used for empirical and simulated data
data <- readRDS('data.RDS')
inits <- list(list(a=rep(0, len=data$n_cross), b=0, sa2=.1, sm2=.05, st2=.03,
              sf2=.15, se2=1, ua=rep(0, len=data$n_animal),
              um=rep(0,len=data$n_ID), ut=rep(0, len=data$n_Tank_ID),
              uf=rep(0, len=data$n_Dam)))
pars <- names(inits[[1]])

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
  verify.models(model=m, params.jags=pars, inits=inits, data=data,
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
        fameffphi_raw=as.numeric(sims.ind[i, grep('fameffphi_raw\\.', x=names(sims.ind))]),
        fameffp_raw=as.numeric(sims.ind[i, grep('fameffp_raw\\.', x=names(sims.ind))]),
        yeareffphi_raw=as.numeric(sims.ind[i, grep('yeareffphi_raw\\.', x=names(sims.ind))])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=pars, inits=inits, data=data,
              lambda=lambda.vec, delta=delta.vec, metric=metric, seeds=seeds,
              Nout=Nout)

message(paste('Finished with model:', m))
setwd('../..')

