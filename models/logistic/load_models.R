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
    logK~dunif(log(1000), log(15000)) #dnorm(log(5000), pow(.21,-.5))
    logr~dunif(log(.01), log(.5)) # dnorm(log(.1), pow(.2,-.5))
    r <- exp(logr)
    K <- exp(logK)
    N[1] <- K
    ## Loop through and project population dynamics
    for(i in 2:num_years){
        temp[i] <- N[i-1]+r*N[i-1]*(1-N[i-1]/K)-catches[i-1]
        N[i] <- max(10/(2-temp[i]/10), temp[i])
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

### ------------------------------------------------------------
## ADMB models

data.admb <- with(data, c(num_years, catches, num_obs, log_pop_obs, years_obs, sd_obs, 999))
setwd('ADMB/')
write.table(data.admb, file='logistic.dat', row.names=FALSE, col.names=FALSE)
system('logistic', ignore.stdout=TRUE)
setwd('..')


### ------------------------------------------------------------
rm(data, init, logistic.traj, data.admb)
message("Finished loading logistic models")

## Verify the model posteriors are the same by getting independent samples
## and comparing.
temp <- jags(data=data.jags, inits=inits.jags, param=params.jags, model.file=model.jags,
     n.chains=1, n.burnin=1000, n.iter=5*10000+1000, n.thin=5)
.jags <- data.frame(temp$BUGSoutput$sims.list)
temp2 <- stan(fit=model.stan, data=data.stan, iter=10100*100, chains=1,
             warmup=1000, thin=100, init=inits.stan,
             control=list(adapt_engaged=TRUE))
.stan <- data.frame(extract(temp2, permuted=FALSE)[,1,])[,c('r', 'K', 'lp__')]
tune.pars <- data.frame(get_sampler_params(temp2))
tune.pars$iteration <- 1:nrow(tune.pars)
ggplot(melt(tune.pars, 'iteration'), aes(iteration, value))+
    facet_wrap('variable', scales='free') + geom_point()
setwd('ADMB/')
file.remove('logistic.psv')
system('logistic -mcmc 1e7 -mcscale 1000 -nosdmcmc -mcsave 10000 -mcseed 5', ignore.stdout=TRUE)
temp3 <- exp(R2admb::read_psv('logistic'))[-1,]
.admb <- data.frame(r=temp3[,1], K=temp3[,2])
setwd('..')


par(mfrow=c(3,3))
acf(.stan$r)
acf(.stan$K)
plot(.stan$r, .stan$K, xlim=c(0,.5), ylim=c(4000, 15000))
acf(.jags$r)
acf(.jags$K)
plot(.jags$r, .jags$K, xlim=c(0,.5), ylim=c(4000, 15000))
qqplot(.stan$r, .jags$r); abline(0,1)
qqplot(.stan$K, .jags$K); abline(0,1)
qqplot(-.stan$lp__, .jags$deviance); abline(0,1)

par(mfrow=c(3,3))
acf(.admb$r)
acf(.admb$K)
plot(.admb$r, .admb$K, xlim=c(0,.5), ylim=c(4000, 15000))
acf(.jags$r)
acf(.jags$K)
plot(.jags$r, .jags$K, xlim=c(0,.5), ylim=c(4000, 15000))
qqplot(.admb$r, .jags$r); abline(0,1)
qqplot(.admb$K, .jags$K); abline(0,1)
qqplot(-.admb$lp__, .jags$deviance); abline(0,1)
