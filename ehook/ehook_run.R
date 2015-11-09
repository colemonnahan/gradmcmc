## An analysis of the Hamley and Skud (1978) data on the effect of hook
## spacing.

## The original data, to be processed into the correct format for each
## model below
hs <- read.csv('data/hs_data.csv')
spacing <- hs$spacing
group <- as.numeric(hs$group)
## The dependent variable (catch.per.hook for now) on log scale since it
## needs to be positive
yobs <- hs$catch.per.hook
log.yobs <- log(yobs+.1)
Nobs <- length(yobs)
Ngroup <- length(unique(group))
day <- hs$daynumber

## Got some reasonable values from an initial run
inits <- list(beta=.1, gamma=.05, logcpue=rep(1,len=Ngroup),
              logcpue_mean=.2, logcpue_sd=.5, logsigma_obs=rep(-.5, len=Ngroup),
              sigma_obs_mean=-.5, sigma_obs_sd=.3)


## These are global parameters. Each software platform specifies them
## slightly different so see the function calls.
n.out1 <- 100
n.thin1 <- 2
n.chains1 <- 1
n.iter1 <- 1.25*n.out1*n.thin1
n.burnin1 <- .2*n.iter1

n.out.ind <- 1000
n.thin.ind <- 1000
n.chains.ind <- 4
n.iter.ind <- 1.25*n.out.ind*n.thin.ind
n.burnin.ind <- .2*n.iter.ind

### ------------------------------------------------------------
## Run the JAGS models
data.jags <-
    list(log_yobs=log.yobs, group=group, Nobs=Nobs, Ngroup=Ngroup, day=day,
         spacing=spacing)
params.jags <-
    c("logcpue_mean", "logcpue_sd", "logsigma_obs", "gamma", "logcpue",
      "sigma_obs_mean", "sigma_obs_sd", "beta")
## JAGS models
model.jags <- function(){
### The Hamley and Skud formula with a random effect on CPUE. Here there is
### a random effect on the observation error
    ## hyperpriors
    logcpue_mean~dunif(-5,5)            # Mean CPUE (asymptote)
    logcpue_sd~dunif(0,5)             # SD of mean CPUE
    sigma_obs_mean~dunif(-1,5)               # mean Observation error, log scale
    sigma_obs_sd~dunif(0,3)               # sd of obs errors, log scale
    beta~dunif(0,10)
    ## priors
    gamma~dunif(0,1)                   # Impact of day on CPUE
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Ngroup){
        logcpue[i]~dnorm(logcpue_mean, pow(logcpue_sd, -2))
        logsigma_obs[i]~dnorm(sigma_obs_mean, pow(sigma_obs_sd, -2))
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        ypred[i] <- exp(logcpue[group[i]])*exp(-day[i]*gamma)*(1-exp(-beta*spacing[i]))
        ## Likelihood of data
        log_yobs[i]~dnorm(log(ypred[i]), pow(exp(logsigma_obs[group[i]]), -2))
    }
}

## Run a long chain with thinning to get independent samples to make sure
## the models are matching.
results.jags.independent <- jags(data=data.jags, parameters.to.save=params.jags,
                           model.file=model.jags, n.chains=n.chains.ind, n.iter=n.iter.ind,
                           n.burnin=n.burnin.ind, n.thin=n.thin.ind)
saveRDS(results.jags.independent, file='results/results.jags.independent.RDS')

results.jags <- jags(data=data.jags, parameters.to.save=params.to.save.jags,
                           model.file=mod.jags, n.chains=n.chains1, n.iter=n.iter1,
                           n.burnin=n.burnin1, n.thin=n.thin1)
saveRDS(results1.jags, file='results/results1.jags.RDS')


## End of JAGS runs for this model
### ------------------------------------------------------------

### ------------------------------------------------------------
## Stan runs
data.stan <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
##if(!exists('stan.model'))
stan.model <- stan(file='ehook.stan', data=data.stan, iter=1, chains=1, init=list(inits))
results1.stanHMC1 <-
    stan(fit=stan.model, data=data.stan, iter=1000, warmup=50,
         chains=n.chains1, thin=1, algorithm='HMC', init=list(inits),
         seed=11212, control=list(adapt_engaged=TRUE, int_time=15))

results2.stan <-
    stan(fit=stan.model, data=data.stan, iter=1000, warmup=50,
         chains=n.chains2, thin=1, algorithm='HMC', init=list(inits),
         seed=11212, control=list(adapt_engaged=TRUE, int_time=15))
saveRDS(results1.stan, file='results/results1.stan.RDS')


### ------------------------------------------------------------
## devtools::install_github('kaskr/adcomp/TMB')
library(TMB)
## TMB runs
data.tmb <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
data.tmb$group <- data.tmb$group-1
## Need to massage inits b/c of transformed parameters
inits.tmb <- inits
boundpinv <- function(x, min, max){
    -log( (max-min)/(x-min) -1)
}
inits.tmb$logcpue_mean2=boundpinv(inits.tmb$logcpue_mean, -5, 5)
inits.tmb$logcpue_sd2=boundpinv(inits.tmb$logcpue_sd, 0, 5)
inits.tmb$sigma_obs_mean2=boundpinv(inits.tmb$sigma_obs_mean, -5, 5)
inits.tmb$sigma_obs_sd2=boundpinv(inits.tmb$sigma_obs_sd, 0,5)
inits.tmb$beta2= boundpinv(inits.tmb$beta, 0,10)
inits.tmb$gamma2=boundpinv(inits.tmb$gamma, 0,1)
inits.tmb$logcpue_mean=NULL
inits.tmb$logcpue_sd=NULL
inits.tmb$sigma_obs_mean=NULL
inits.tmb$sigma_obs_sd=NULL
inits.tmb$beta=NULL
inits.tmb$gamma=NULL

## Need to reorder parameter list for TMB.

pars <- c("beta2","gamma2","logcpue", "logcpue_mean2","logcpue_sd2",
          "logsigma_obs", "sigma_obs_mean2", "sigma_obs_sd2" )
inits.tmb <- inits.tmb[pars]
compile("ehook.cpp")
dyn.load(TMB::dynlib("ehook"))
tmb.model <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='ehook')
tmb.model.opt <- do.call(optim, tmb.model)

results.tmb.independent <-
    TMB::mcmc(obj=tmb.model, nsim=1500, algorithm='NUTS', Madapt=500 )
results.tmb.independent$LP <- -apply(results.tmb.independent, 1, tmb.model$fn)
saveRDS(results.tmb.independent, file='results/results.tmb.independent.RDS')
