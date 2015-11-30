
data.logistic <-
    logistic.traj(logr=logr.true, logK=logK.true, num.years=100,
                  sd.catch=.5, prop.caught=.05, years.obs=c(15,35,85,90),
                  sd.obs=.1)
data.logistic.TMB <- data.logistic

data.jags <- data.logistic
data.jags$log_pop_obs <- log(data.jags$pop_obs)
data.jags$pop_obs <- NULL
params.jags <- c("r", "K", "N")
model.jags <- function(){
    ## Priors
    logK~dnorm(log(5000), pow(.21,-2))
    logr~dnorm(log(.1), pow(.2,-2))
    r <- exp(logr)
    K <- exp(logK)
    N[1] <- K
    ## Loop through and project population dynamics
    for(i in 2:num_years){
        N[i] <- max(10,N[i-1]+r*N[i-1]*(1-N[i-1]/K)-catches[i])
    }
    ## Likelihood of data
    for(j in 1:num_obs){
        log_pop_obs[j]~dnorm(log(N[years_obs[j]]), pow(sd_obs, -2))
        ## log_pop_obs[j]~dnorm(8.5, pow(sd_obs, -2))
    }
}

temp <- jags(data=data.jags, inits=list(list(logr=log(.2), logK=log(5000))), param=params.jags, model.file=model.jags,
     n.chains=1, n.burnin=1000, n.iter=10050, n.thin=1)
xx <- data.frame(temp$BUGSoutput$sims.list)
plot(xx$r, xx$K)
