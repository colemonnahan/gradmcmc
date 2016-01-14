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

## First get independent samples to verify the models are the same
verify.models(params.jag=c('mu'), model.jag='mvn.jags',
              inits.jags=inits.jags, data.jags=data.jags,
              model.stan='mvn.stan', inits.stan=inits.stan,
              data.stan=data.stan, Niter=1e6, Nthin=1000)

## Now rerun across gradient of acceptance rates and compare to JAGS
model.stan <- stan(file='growth.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)

run.chains(model=m, seeds=1, Nout=500, L=NULL,
