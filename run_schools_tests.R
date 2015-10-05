model <- 'schools'
setwd(model)
## playing with simple STAN models
library(rstan)
## this is the default example
schools_dat <- list(J=8, y=c(28,8,-3,7,-1,1,18,12),
                      sigma=c(15,10,16,11,9,11,10,18))
## with(schools_dat, plot(y, sigma))
fit <- stan(file='schools.stan', data=schools_dat, iter=50000, chains=1)
pairs(fit, pars=c('mu', 'tau', 'lp__'))
print(fit)
la <- extract(fit, permuted=TRUE)
mu.stan <- la$mu
tau.stan <- la$tau
mu <- 7.87
tau <- 6.67
logtau <- log(tau)
theta <- c(11.39,8,6,8,5,7,10,8)
eta <- (theta-mu)/exp(logtau)


## Now fit it in JAGS too
library(R2jags)
data.jags <- list(y=schools_dat$y, sigma=schools_dat$sigma)
params.to.save <- c("mu", "tau", "eta")
## JAGS model
model <- function(){  ## The 8-schools problem in JAGS
    ## priors for fixed effects
    mu~dunif(-20,20)
    tau~dunif(0,100)
    ## Loop through observations and calculate likelihood
    for(i in 1:8){
        ## Hyperdistribution for eta
        eta[i]~dnorm(0,1)
        ## The predicted values
        theta[i] <- mu+ tau*eta[i]
        ## Likelihood of data
        y[i]~dnorm(theta[i], pow(sigma[i],-2))
    }
}
jags.fit <- jags(data=data.jags, inits=list(list(mu=mu, tau=exp(logtau), eta=eta)),
             parameters.to.save=params.to.save,
             model.file=model,
             n.chains=1, n.iter=10500, n.burnin=500,
             n.thin=1)
jags.pars <- as.data.frame(jags.fit$BUGSoutput$sims.matrix)
## par(mfrow=c(3,4))
## for(i in 1:ncol(jags.pars)) acf(jags.pars[,i])
apply(jags.pars, 2, mean)
tau.jags <-  jags.pars$tau
qqplot(tau.stan, tau.jags); abline(0,1)


## recreate 8-schools in TMB
## -sum(dnorm(x=theta, mean=schools_dat$y, sd=schools_dat$sigma,
##                   log=TRUE)) - sum(dnorm(eta, 0,1, log=TRUE))
## init <- function(){list(mu=mu, logtau=logtau, eta=eta)}
## stan(file = 'schools.stan', data = schools_dat,
##             iter = 1, chains = 1, init=init)
library(TMB)
compile("schools.cpp")
dyn.load(TMB::dynlib("schools"))
data <- list(Y=schools_dat$y, sigma=schools_dat$sigma)
parameters <- list(mu=mu, logtau=2, eta=eta)
## ## Use an iteration from JAGS to see if LP matches.
## parameters <- list(mu=jags.pars[5,'mu'], logtau=log(jags.pars[5,'tau']),
##                    eta=as.numeric(jags.pars[5, 2:9]))
schools <- TMB::MakeADFun(data,parameters, DLL='schools')
schools$env$beSilent()
temp <- llply(xseq <- seq(-5,5, len=100), function(x) schools$fn(c(mu, x, eta)))
plot(50/(1+exp(-xseq)), temp)

## get a covariance from JAGS
temp <- jags.pars[, c(10,11, 2:9)]; temp$tau <- log(temp$tau)
schools.covar <- cov(temp)
## schools.mcmc <- run_mcmc(schools, nsim=110, algorithm='NUTS', max_doublings=3,
##                          diagnostic=FALSE, delta=.5, Madapt=100)
## schools.mcmc <- schools.mcmc[-(1:1000),]
schools.mcmc <- run_mcmc(schools, nsim=50000, algorithm='HMC', L=5, eps=.25,
                         diagnostic=FALSE, covar=schools.covar)
## schools.mcmc <- run_mcmc(schools, nsim=10000, algorithm='RWM',
##                          diagnostic=FALSE, alpha=2, covar=schools.covar)
## schools.mcmc$tau <- 50/(1+exp(-schools.mcmc$tau_unbounded))
schools.mcmc$tau <- exp(schools.mcmc$logtau)
schools.mcmc$logtau <- NULL
pairs(schools.mcmc[, c('mu', 'tau')], pch=16, col=rgb(0,0,0,.01))
tau.tmb <- schools.mcmc$tau
mu.tmb <- schools.mcmc$mu
acf(tau.jags)
acf(tau.tmb)
acf(tau.stan)
plot(tau.tmb)
qqplot(tau.stan, tau.tmb);abline(0,1)
