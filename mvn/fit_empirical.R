### ------------------------------------------------------------
### Fit the empirical data

## Load data
Npar <- 5
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))

## Inits
inits <- list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2)
inits.jags <- inits.stan <- list(inits)
data.jags <- data.stan <- data
params.jags <- 'mu'

## First get independent samples to verify the models are the same
verify.models(params.jag=params.jags, model.jag='mvn.jags',
              inits.jags=inits.jags, data.jags=data.jags,
              model.stan='mvn.stan', inits.stan=inits.stan,
              data.stan=data.stan, Niter=1e6, Nthin=1000)

## model.jags <- jags(model.file='mvn.jags', inits=inits.jags, para=params.jags,
##                    n.chains=1, data=data.jags)
## model.stan <- stan(file='mvn.stan', data=data.stan, iter=50, chains=1, alg='HMC',
##                    warmup=10, thin=1, init=inits.stan, control=list(metric='unit_e'))
## xx <- as.data.frame(get_sampler_params(model.stan))

## Now rerun across gradient of acceptance rates and compare to JAGS
results.empirical <-
    run.chains(model=m, seeds=1:3, Nout=10000, lambda=c(.1,1),
               metric=c('diag_e', 'dense_e'),
               delta=c(.3,.5,.8, .95), data.jags=data.jags,
               Nthin=1, inits.jags=inits.jags, params.jags=params.jags,
               data.stan=data.stan, inits.stan=inits.stan)
perf <- results.empirical$perf
adapt <- results.empirical$adapt
plot.empirical.results(perf, adapt)
