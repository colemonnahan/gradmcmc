## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
data <- readRDS('data.RDS')
params.jags <-
    c("yearInterceptSD", "plantInterceptSD", "plantSlopeSD", "intercept",
      "slope", "yearInterceptEffect", "plantSlopeEffect", "plantInterceptEffect")
inits <-
    list(list(
        yearInterceptSD = 1,
              plantInterceptSD = 1,
              plantSlopeSD = 1,
              intercept = rep(0,data$Nstage), slope = 0,
              yearInterceptEffect = rep(0, data$Nyear),
              plantInterceptEffect = rep(0, data$Nplant),
              plantSlopeEffect = rep(0, data$Nplant)))
## stan.fit <- stan(file='wildflower.stan', data=data, init=inits,seed=11,
##                  pars=params.jags, iter=2000, chains=1)
## shinystan::launch_shinystan(stan.fit)

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

message(paste('Finished with model:', m))
setwd('../..')

