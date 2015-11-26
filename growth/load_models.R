message("Loading growth models into the workspace")

## Global parameters
Linf.mean <- 50
k.mean <- .1
t0 <- 0
Linf.sigma <- .01
k.sigma <- .01
sigma.obs <- .1
Ntime <- 50

## source("generate_data.R")
dat <- readRDS('growth_data_50000.RDS')

dat2 <- dat[dat$fish<=Nfish,]
data <- list(Nfish=Nfish, Nobs=nrow(dat2), loglengths=log(dat2$lengths),
                  fish=dat2$fish, ages=dat2$ages)
init <- list(Linf_mean=Linf.mean, Linf_sigma=Linf.sigma,
                  k_mean=k.mean, k_sigma=k.sigma, sigma_obs=sigma.obs,
                  logLinf=log(rep(Linf.mean, len=Nfish)),
                  logk=log(rep(k.mean, len=Nfish)))

### ------------------------------------------------------------
## JAGS models
data.jags <- data
params.jags <-
    c("Linf_mean", "Linf_sigma", "k_mean", "k_sigma", "sigma_obs", "logLinf",
      "logk")
inits.jags <- list(init)
model.jags <- function(){
    ## hyperpriors
    Linf_mean~dunif(30,70)
    Linf_sigma~dunif(0,.05)
    k_mean~dunif(0,.5)
    k_sigma~dunif(0,.05)
    ## priors
    sigma_obs~dunif(0,1)
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Nfish){
        logLinf[i]~dnorm(log(Linf_mean), pow(Linf_sigma, -2))
        logk[i]~dnorm(log(k_mean), pow(k_sigma, -2))
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        ypred[i] <- logLinf[fish[i]]+
            log( (1-exp(-exp(logk[fish[i]])*ages[i])))
        ## Likelihood of data
        loglengths[i]~dnorm(ypred[i], pow(sigma_obs, -2))
    }
}

## temp <- jags(data=data.jags, inits=list(init), param=params.jags, model.file=model.jags,
##      n.chains=1, n.burnin=1000, n.iter=5000, n.thin=5)
## xx <- data.frame(temp$BUGSoutput$sims.list)
## acf(xx$Linf_mean)
## mean(xx$Linf_mean)
## acf(xx$k_mean)
## mean(xx$k_mean)
## acf(xx$Linf_sigma)
## mean(xx$Linf_sigma)
## acf(xx$k_sigma)
## mean(xx$k_sigma)
## pairs(xx[, 1:10])
## End of JAGS
### ------------------------------------------------------------



### ------------------------------------------------------------
## Stan models
data.stan <- data
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
inits.stan <- list(init)
model.stan <- stan(file='growth.stan', data=data.stan, iter=10, chains=1,
                   init=inits.stan)

## End of Stan
### ------------------------------------------------------------

message("Finished loading growth models")


## model.R <- function(Linf_mean, Linf_sigma, k_mean, k_sigma, sigma_obs,
##                     Linf, k){
##     ages <- data$ages; fish <- data$fish; Nobs=data$Nobs
##     loglengths <- data$loglengths
##     NLL <- 0
##     ## hyperpriors
##     NLL <- NLL-dunif(x=Linf_mean, 30,70, log=TRUE)
##     NLL <- NLL-dunif(x=Linf_sigma, 0,0.05, log=TRUE)
##     NLL <- NLL-dunif(x=k_mean, 0, .5, log=TRUE)
##     NLL <- NLL-dunif(x=k_sigma, 0,0.05, log=TRUE)
##     ## priors
##     NLL <- NLL-dunif(sigma_obs, 0,1)
##     ## Loop through the hyperparameters (on group) and calculate
##     ## probabilities
##     NLL <- NLL- sum(dnorm(Linf, Linf_mean, Linf_sigma, log=TRUE))
##     NLL <- NLL- sum(dnorm(k, k_mean, k_sigma, log=TRUE))
##     ypred <- rep(NA, len=Nobs)
##     for(i in 1:Nobs){
##         ypred[i] <- log(Linf[fish[i]])+log( (1-exp(-k[fish[i]]*ages[i])) )
##         ## Likelihood of data
##         NLL <- NLL-dnorm(loglengths[i], 1.5, sigma_obs, log=TRUE)
##     }
##     NLL
## }
## model.R(Linf.mean, Linf.sigma, k.mean, k.sigma, sigma.obs, Linf, k)
