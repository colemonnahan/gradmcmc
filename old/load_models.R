### ------------------------------------------------------------
## The purpose of this file is to prepare the models for running later. May
## need to refresh/create data or not depending on the model. The models
## are loaded into the workspace and used later on.

## Generate data
ss_logistic.traj <- function(r, K, num.years, sd.catch, prop.caught,
                             years.obs, sd.obs, sd.process, F, plot=TRUE){
    catches <- trajectory <- rep(0,num.years)
    u <- rnorm(n=num.years, mean=0, sd=sd.process)
    trajectory[1] <- K + u[1]
    for( yr in 2: num.years){
        Nt <- trajectory[yr-1]
        Ct <- catches[yr-1] <- Nt*(1-exp(-ifelse(Nt<K/10, 0, F)))
        trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct) * exp(u[yr]-sd.process^2/2)
        if( trajectory[yr]<=0 ) { trajectory[yr] <- NA; break}
    }
    ## generate lognormal samples
    log.pop.obs <-
        rnorm(n=length(years.obs),
              mean=log(trajectory[years.obs]),sd=sd.obs) -sd.obs^2/2
    if(plot){
        plot(1:num.years, ylim=c(0,K*1.1),y= trajectory, type='l')
        points(1:num.years, catches, type='h')
        points(years.obs, exp(log.pop.obs), col='red')
    }
    return(list(N=num.years, catches=catches, u=u, logcpue=log.pop.obs))
    ## return(trajectory)
}
r.true <- .1
K.true <- 5000
sd.process <- .02
sd.obs <- .02
data <- ss_logistic.traj(r=r.true, K=K.true, num.years=Nyears,
                         years.obs=1:Nyears, sd.obs=sd.obs,
                         sd.process=sd.process, F=.05,
                         plot=FALSE)
init <- list(r=r.true, K=K.true, q=1, sd_obs=.1, sd_process=sd.process,
             u=data$u)
u <- data$u
data$u <- NULL
## check initial values
r <- (init$r)
K <- (init$K)
trajectory <- rep(NA, data$N)
trajectory[1] <- K*exp(init$u[1]-sd.process^2/2)
for( yr in 2: data$N){
    Nt <- trajectory[yr-1]
    Ct <- data$catches[yr-1]
    trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct)*exp(init$u[yr]-sd.process^2/2)
}
plot(trajectory, ylim=c(0,max(exp(data$logcpue))), type='l')
with(data, points(exp(logcpue), col='red'))

### ------------------------------------------------------------
## JAGS models
data.jags <- data
params.jags <- c('r', 'K', 'sd_obs', 'q', 'sd_process', 'u')
inits.jags <- list(init)
model.jags <- 'ss_logistic.jags'
## temp <- jags(data=data.jags, inits=list(init), param=params.jags, model.file='ss_logistic.jags',
##             n.chains=1, n.burnin=0, n.iter=5, n.thin=1)
## End of JAGS
### ------------------------------------------------------------

### ------------------------------------------------------------
## Stan models
data.stan <- data
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
inits.stan <- list(init)
model.stan <- stan(file='ss_logistic.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)

## rm(data, init,  K.true, r.true, sd.process)
## End of Stan
### ------------------------------------------------------------
message("Finished loading ss_logistic models")

## Some development code to test the models.
thin <- 1000
warmup <- 5000
Niter <- 1000
pars <- c('r', 'K', 'q', 'sd_obs', 'sd_process')
temp <- jags(data=data.jags, inits=inits.jags, param=params.jags, model.file='ss_logistic.jags',
     n.chains=1, n.burnin=warmup, n.iter=thin*Niter+warmup, n.thin=thin)
.jags <- data.frame(temp$BUGSoutput$sims.list)
par(mfrow=c(3,5), mar=c(2,2,2,0))
trash <- lapply(pars, function(x)
    hist(.jags[,x], breaks=30, main=x))
trash <- lapply(pars, function(x)
    acf(.jags[,x], main=x))
trash <- lapply(pars, function(x)
    plot(.jags[,x], type='l',main=x))
temp2 <- stan(fit=model.stan, data=data.stan, iter=Niter*thin+warmup, chains=1,
             warmup=warmup, thin=thin, init=inits.stan,
             control=list(adapt_engaged=TRUE, max_treedepth=7))

.stan <- data.frame(extract(temp2, permuted=FALSE)[,1,])
tune.pars <- data.frame(get_sampler_params(temp2))
tune.pars$iteration <- 1:nrow(tune.pars)
tune.pars$stepsize__ <- log(tune.pars$stepsize__)
ggplot(melt(tune.pars, 'iteration'), aes(iteration, value))+
    facet_wrap('variable', scales='free') + geom_point()

par(mfrow=c(3,5), mar=c(2,2,2,0))
trash <- lapply(pars, function(x)
    hist(.stan[,x], breaks=30, main=x))
trash <- lapply(pars, function(x)
    acf(.stan[,x], main=x))
trash <- lapply(pars, function(x)
    plot(.stan[,x], type='l',main=x))

pairs(.stan[,pars])
pairs(.jags[,pars])

par(mfrow=c(3,6))
acf(.stan$r)
acf(.stan$K)
acf(.stan$q)
acf(.stan$sd_process)
acf(.stan$sd_obs)
plot(.stan$r, .stan$K, xlim=c(0,.2), ylim=c(2000, 10000))
acf(.jags$r)
acf(.jags$K)
acf(.jags$q)
acf(.jags$sd_process)
acf(.jags$sd_obs)
plot(.jags$r, .jags$K, xlim=c(0,.2), ylim=c(2000, 10000))
qqplot(.stan$r, .jags$r); abline(0,1)
qqplot(.stan$K, .jags$K); abline(0,1)
qqplot(.stan$q, .jags$q); abline(0,1)
qqplot(.stan$sd_process, .jags$sd_process); abline(0,1)
qqplot(.stan$sd_obs, .jags$sd_obs); abline(0,1)
qqplot(-.stan$lp__, .jags$deviance); abline(0,1)

jags.ess <- effectiveSize(.jags)
stan.ess <- effectiveSize(.stan)
plot(jags.ess, stan.ess)

par(mfrow=c(5,10), mar=c(0,0,0,0))
for(i in 1:Nyears){
    plot(.stan[, paste0('u.',i,'.')], type='l')
    abline(h=u[i], col='red')
}
par(mfrow=c(5,10), mar=c(0,0,0,0))
for(i in 1:Nyears){
    plot(.jags[, paste0('u.',i)], type='l')
    abline(h=u[i], col='red')
}
par(mfrow=c(5,10), mar=c(0,0,0,0))
for(i in 1:Nyears){
    qqplot(.jags[, paste0('u.',i)], .stan[, paste0('u.',i,'.')],
           axes=FALSE)
    box()
    abline(0,1, col='red')
}

## Plot model fits for both softwares
get.trajectory <- function(data, par){
    u <- as.numeric(par[,grep('u.', x=names(par))])
    trajectory <- rep(NA, data$N)
    trajectory[1] <- par$K*exp(u[1])
    for(yr in 2:data$N){
        Nt <- trajectory[yr-1]
        Ct <- data$catches[yr-1]
        trajectory[yr] <- (Nt+par$r*Nt*(1- (Nt/par$K))-Ct)*exp(u[yr])
    }
    return(trajectory)
}

## check model fits
r <- (init$r)
K <- (init$K)
trajectory <- rep(NA, data$N)
trajectory[1] <- K*exp(init$u[1]-sd.process^2/2)
for( yr in 2: data$N){
    Nt <- trajectory[yr-1]
    Ct <- data$catches[yr-1]
    trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct)*exp(init$u[yr]-sd.process^2/2)
}
plot(trajectory, ylim=c(0,max(exp(data$logcpue))), type='l', lwd=2)
for(i in 1:100){
    lines(1:data$N, get.trajectory(data=data, par=.jags[i,]), col=gray(.4))
}
lines(trajectory, lwd=2)
with(data, points(exp(logcpue), col='red'))

## Check distribution of biomass in each year
.jags.temp <- ldply(1:nrow(.jags), function(i)
    get.trajectory(data=data, par=.jags[i,]))
.stan.temp <- ldply(1:nrow(.stan), function(i)
    get.trajectory(data=data, par=.stan[i,]))
par(mfrow=c(5,10), mar=c(0,0,0,0))
for(i in 1:data$N){
    qqplot(.jags.temp[,i], .stan.temp[,i], axes=FALSE)
    box()
    abline(0,1, col='red')
}
