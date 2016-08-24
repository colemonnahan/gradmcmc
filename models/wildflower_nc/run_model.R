## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
data <- readRDS('data.RDS')
params.jags <-
    c("yearInterceptSD", "plantInterceptSD", "plantSlopeSD", "intercept",
      "slope", "yearInterceptEffect_raw", "plantSlopeEffect_raw", "plantInterceptEffect_raw")
inits <-
    list(list(
        yearInterceptSD = 1,
              plantInterceptSD = 1,
              plantSlopeSD = 1,
              intercept = rep(0,data$Nstage), slope = 0,
              yearInterceptEffect_raw= rep(0, data$Nyear),
              plantInterceptEffect_raw= rep(0, data$Nplant),
              plantSlopeEffect_raw= rep(0, data$Nplant)))

## stan.fit <- stan(file='wildflower_nc.stan', data=data, init=inits,seed=11,
##                  pars=params.jags, iter=1000, chains=1)
## shinystan::launch_shinystan(stan.fit)
## jags.fit <- jags(data=data, inits=inits, parameters.to.save=params.jags,
##                  model.file='wildflower_nc.jags', n.chains=1, n.iter=2000)
## jags.sims <- shinystan::as.shinystan(jags.sims)
## shinystan::launch_shinystan(jags.sims)

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
  verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
              Nout=Nout.ind, Nthin=Nthin.ind)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]

inits <- lapply(1:length(seeds), function(i)
    list(
        yearInterceptSD=sims.ind$yearInterceptSD[i],
        plantInterceptSD=sims.ind$plantInterceptSD[i],
        plantSlopeSD=sims.ind$plantSlopeSD[i],
        intercept = as.numeric(sims.ind[i, grep('intercept', names(sims.ind))]),
        yearInterceptEffect =
            as.numeric(sims.ind[i, grep('yearInterceptEffect', names(sims.ind))]),
        plantInterceptEffect =
            as.numeric(sims.ind[i, grep('plantInterceptEffect', names(sims.ind))]),
        plantSlopeEffect =
            as.numeric(sims.ind[i, grep('plantSlopeEffect', names(sims.ind))]),
        slope=sims.ind$slope[i])
)

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

## library(coda)
## library(shinystan)
## stan.fit <- readRDS(file='fits/stan_nuts_diag_e_0.8_10_.RDS')
## stan.fit <- as.shinystan(mcmc.list(as.mcmc(readRDS(file='sims.ind.RDS'))))
## shinystan::launch_shinystan(stan.fit)


message(paste('Finished with model:', m))
setwd('../..')

