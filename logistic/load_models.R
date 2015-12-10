set.seed(115)
message('Generating data for logistic model')
## Generate fake data
logistic.traj <- function(logr, logK, num.years, sd.catch, prop.caught,
                          years.obs, sd.obs){
    r <- exp(logr)
    K <- exp(logK)
    catches <- trajectory <- rep(0,num.years)
    trajectory[1] <- K
    for( yr in 2: num.years){
        Nt <- trajectory[yr-1]
        Ct <- catches[yr-1] <-
            if(yr<50) exp(rnorm(1,log(prop.caught*Nt),sd=sd.catch)) else 0
        trajectory[yr] <- Nt+r*Nt*(1- (Nt/K))-Ct
        if( trajectory[yr]<=0 ) { trajectory[yr] <- NA; break}
    }
    ## generate lognormal samples
    log.pop.obs <-
        rnorm(n=length(years.obs),
                  mean=log(trajectory[years.obs]),sd=sd.obs) -sd.obs^2/2
    plot(1:num.years, ylim=c(0,K*1.1),y= trajectory, type='b')
    points(1:num.years, catches, type='h')
    points(years.obs, exp(log.pop.obs), col='red')
    return(list(num_years=num.years, catches=catches,
                num_obs=length(log.pop.obs), log_pop_obs=log.pop.obs,
                sd_obs=sd.obs, years_obs=years.obs))
    ## return(trajectory)
}
logr.true <- log(.1)
logK.true <- log(5000)
data <-
    logistic.traj(logr=logr.true, logK=logK.true, num.years=100,
                  sd.catch=.5, prop.caught=.05, years.obs=c(65, 75, 85,90),
                  sd.obs=.05)
init <- list(logr=logr.true, logK=logK.true)
### ------------------------------------------------------------
## JAGS models
data.jags <- data
data.jags$pop_obs <- NULL
params.jags <- c("r", "K")
inits.jags <- list(init)
model.jags <- function(){
    ## Priors
    logK~dunif(log(2000), log(15000)) #dnorm(log(5000), pow(.21,-.5))
    logr~dunif(log(.01), log(.5)) # dnorm(log(.1), pow(.2,-.5))
    r <- exp(logr)
    K <- exp(logK)
    N[1] <- K
    ## Loop through and project population dynamics
    for(i in 2:num_years){
        N[i] <- max(1,N[i-1]+r*N[i-1]*(1-N[i-1]/K)-catches[i])
    }
    ## Likelihood of data
    for(j in 1:num_obs){
        log_pop_obs[j]~dnorm(log(N[years_obs[j]]), pow(sd_obs, -2))
        ## log_pop_obs[j]~dnorm(8.5, pow(sd_obs, -2))
    }
}


### ------------------------------------------------------------
## Stan models
data.stan <- data
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
inits.stan <- list(init)
model.stan <- stan(file='logistic.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)
## End of Stan
### ------------------------------------------------------------
rm(data, init, logistic.traj)
message("Finished loading logistic models")

temp <- jags(data=data.jags, inits=inits.jags, param=params.jags, model.file=model.jags,
     n.chains=1, n.burnin=1000, n.iter=2000, n.thin=1)
xx <- data.frame(temp$BUGSoutput$sims.list)
temp2 <- stan(fit=model.stan, data=data.stan, iter=10*1000, chains=1,
             warmup=1000, thin=10, init=inits.stan,
             control=list(adapt_engaged=TRUE))
yy <- data.frame(extract(temp2, permuted=FALSE)[,1,])[,c('r', 'K')]
pars <- data.frame(get_sampler_params(temp))
pars$iteration <- 1:nrow(pars)
ggplot(melt(pars, 'iteration'), aes(iteration, value))+
    facet_wrap('variable', scales='free') + geom_point()
par(mfrow=c(2,3))
acf(yy$r)
acf(yy$K)
plot(yy$r, yy$K)
acf(xx$r)
acf(xx$K)
plot(xx$r, xx$K)

## Check the posteriors are the same
par(mfrow=c(2,3))
qqplot(yy$r, xx$r); abline(0,1)
qqplot(yy$K, xx$K); abline(0,1)
