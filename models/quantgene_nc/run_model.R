## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
## The true values used for empirical and simulated data
data <- readRDS('data.RDS')
inits <- list(list(a=rep(0, len=data$n_cross), b=0, sa2=.1, sm2=.05, st2=.03,
              sf2=.15, se2=1, ua_raw=rep(0, len=data$n_animal),
              um_raw=rep(0,len=data$n_ID), ut_raw=rep(0, len=data$n_Tank_ID),
              uf_raw=rep(0, len=data$n_Dam)))
pars <- names(inits[[1]])

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
  verify.models(model=m, params.jags=pars, inits=inits, data=data,
                Nout=Nout.ind, Nthin=Nthin.ind)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i)
    list(a=as.numeric(sims.ind[i, grep('a\\.', x=names(sims.ind))]),
         b=sims.ind$b[i], sa2=sims.ind$sa2[i], sm2=sims.ind$sm2[i],
         st2=sims.ind$st2[i], sf2=sims.ind$sf2[i], se2=sims.ind$se2[i],
      ua_raw=as.numeric(sims.ind[i, grep('ua_raw\\.', x=names(sims.ind))]),
      um_raw=as.numeric(sims.ind[i, grep('um_raw\\.', x=names(sims.ind))]),
      ut_raw=as.numeric(sims.ind[i, grep('ut_raw\\.', x=names(sims.ind))]),
      uf_raw=as.numeric(sims.ind[i, grep('uf_raw\\.', x=names(sims.ind))])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=pars, inits=inits, data=data,
              lambda=lambda.vec, delta=c(.9) , metric=metric,
              seeds=seeds, Nout=Nout)

message(paste('Finished with model:', m))
setwd('../..')

