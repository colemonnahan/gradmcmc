message("Loading ehook models into the workspace")

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

### ------------------------------------------------------------
## TMB models
data.tmb <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
data.tmb$group <- data.tmb$group-1
## Need to massage inits b/c of transformed parameters
inits.tmb <- inits
inits.tmb$logcpue_mean2=boundpinv(inits.tmb$logcpue_mean, -5, 5)
inits.tmb$logcpue_sd2=boundpinv(inits.tmb$logcpue_sd, 0, 5)
inits.tmb$sigma_obs_mean2=boundpinv(inits.tmb$sigma_obs_mean, -5, 5)
inits.tmb$sigma_obs_sd2=boundpinv(inits.tmb$sigma_obs_sd, 0,5)
inits.tmb$beta2= boundpinv(inits.tmb$beta, 0,1)
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
model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='ehook')
model.tmb$hessian <- TRUE
model.tmb.opt <- do.call(optim, model.tmb)
covar.tmb <- solve(model.tmb.opt$hessian)
## This is given as a vector but need to put it in list form to put back
## into model
est.tmb <- model.tmb.opt$par
x <- as.vector(est.tmb)
inits.mle <- list(beta2=x[1], gamma2=x[2], logcpue=x[3:16], logcpue_mean2=x[17],
     logcpue_sd2=x[18], logsigma_obs=x[19:32], sigma_obs_mean2=x[33],
     sigma_obs_sd2=x[34])
inits.tmb <- inits.mle

## End of TMB
### ------------------------------------------------------------


### ------------------------------------------------------------
## JAGS models
data.jags <-
    list(log_yobs=log.yobs, group=group, Nobs=Nobs, Ngroup=Ngroup, day=day,
         spacing=spacing)
params.jags <-
    c("logcpue_mean", "logcpue_sd", "logsigma_obs", "gamma", "logcpue",
      "sigma_obs_mean", "sigma_obs_sd", "beta")
model.jags <- function(){
    ## hyperpriors
    logcpue_mean~dunif(-5,5)            # Mean CPUE (asymptote)
    logcpue_sd~dunif(0,5)             # SD of mean CPUE
    sigma_obs_mean~dunif(-1,5)               # mean Observation error, log scale
    sigma_obs_sd~dunif(0,3)               # sd of obs errors, log scale
    beta~dunif(0,1)
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
## End of JAGS
### ------------------------------------------------------------

### ------------------------------------------------------------
## Stan models
data.stan <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
##if(!exists('model.stan'))
inits.stan <- list(inits)
model.stan <- stan(file='ehook.stan', data=data.stan, iter=1, chains=1,
                   init=inits.stan)
## End of Stan
### ------------------------------------------------------------



message("Finished loading ehook models")
