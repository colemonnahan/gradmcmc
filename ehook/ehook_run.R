## An analysis of the Hamley and Skud (1978) data on the effect of hook
## spacing.

## The parameters for the short chain


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

## These are global parameters. Each software platform specifies them
## slightly different so see the function calls.
n.out1 <- 100
n.thin1 <- 2
n.chains1 <- 1
n.iter1 <- 1.25*n.out1*n.thin1
n.burnin1 <- .2*n.iter1

n.out2 <- 1000
n.thin2 <- 100
n.chains2 <- 1
n.iter2 <- 1.25*n.out2*n.thin2
n.burnin2 <- .2*n.iter2

### ------------------------------------------------------------
## Run the JAGS models
data.jags <-
    list(log_yobs=log.yobs, group=group, Nobs=Nobs, Ngroup=Ngroup, day=day,
         spacing=spacing)
params.to.save.jags <-
    c("logcpue_mean", "logcpue_sd", "logsigma_obs", "gamma", "logcpue",
      "sigma_obs_mean", "sigma_obs_sd", "beta")
## JAGS models
mod.jags <- function(){
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
    ## The precision params
    logcpue_tau <- pow(logcpue_sd, -2)
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Ngroup){
        logcpue[i]~dnorm(logcpue_mean, logcpue_tau)
        cpue[i] <- exp(logcpue[i])
        logsigma_obs[i]~dnorm(sigma_obs_mean, pow(sigma_obs_sd, -2))
        sigma_obs[i] <- exp(logsigma_obs[i])
        obs_tau[i] <- pow(sigma_obs[i], -2)
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        ## The predicted values
        ypred[i] <- cpue[group[i]]*exp(-day[i]*gamma)*(1-exp(-beta*spacing[i]))
        ## Likelihood of data
        log_yobs[i]~dnorm(log(ypred[i]), obs_tau[group[i]])
    }
}

## results1.jags.list <- jags(data=data.jags, parameters.to.save=params.to.save.jags,
##                            model.file=mod.jags, n.chains=n.chains1, n.iter=n.iter1,
##                            n.burnin=n.burnin1, n.thin=n.thin1)
## saveRDS(results1.jags.list, file='results/results1.jags.list.RDS')
results2.jags.list <- jags(data=data.jags, parameters.to.save=params.to.save.jags,
                           model.file=mod.jags, n.chains=n.chains2, n.iter=n.iter2,
                           n.burnin=n.burnin2, n.thin=n.thin2)
saveRDS(results2.jags.list, file='results/results2.jags.list.RDS')


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
results1.stan.list <-
    stan(fit=stan.model, data=data.stan, iter=1000, warmup=50,
         chains=n.chains1, thin=1, algorithm='HMC', init=list(inits),
         seed=11212, control=list(adapt_engaged=TRUE, int_time=15))
saveRDS(results1.stan.list, file='results/results1.stan.list.RDS')




