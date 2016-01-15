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

## Now rerun across gradient of acceptance rates and compare to JAGS
model.stan <- stan(file='mvn.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)
xx <- as.data.frame(get_sampler_params(model.stan))

yy <- run.chains(model=m, seeds=2, Nout=1000, L=c(.1,1),
                 delta=c(.3,.4,.5,.8,.9, .95, .99),
                 data.jags=data.jags, Nthin=3,
                 inits.jags=inits.jags, params.jags=params.jags,
                 data.stan=data.stan, inits.stan=inits.stan)

perf2 <- subset(yy$perf, platform!='jags')
perf.jags <- subset(yy$perf, platform=='jags')
ggplot(perf2, aes(delta.final, perf, group=seed))+ geom_line()+
    facet_grid(.~platform) + geom_hline(yintercept=perf.jags$perf)
ggplot(perf2, aes(delta.target, delta.final, color=platform))+
    geom_point() + xlim(0,1) + ylim(0,1)

