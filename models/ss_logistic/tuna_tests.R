### ------------------------------------------------------------
## The purpose of this file is to prepare the models for running later. May
## need to refresh/create data or not depending on the model. The models
## are loaded into the workspace and used later on.
## Generate data
dat <- read.csv('tuna_data.csv')
cpue <- dat$CPUE
catches <- dat$Catches
data <- list(catches=catches, logcpue=log(cpue), N=nrow(dat))
init <- list(logr=log(.8), logK=log(279), iq=5, isigma2=100, itau2=100,
             u=rep(0, len=nrow(dat)))

### ------------------------------------------------------------
## JAGS models
data.jags <- data
params.jags <- c('r', 'K', 'sigma2', 'q', 'tau2', 'B[2]', 'B1990', 'P1990')
inits.jags <- list(init)
model.jags <- function(){
    logK~dnorm(5.04, pow(0.5162,-2))
    logr~dnorm(-1.38, pow(0.51, -2))
    iq~dgamma(0.001, 0.001)
    isigma2~dgamma(3.785518, 0.010223)
    itau2~dgamma(1.708603, 0.008613854)
    K <- exp(logK)
    r <- exp(logr)
    q <- 1/iq
    sigma2 <- pow(isigma2,-1)
    tau2 <- pow(itau2,-1)
    B[1] <- K
    u[1]~dnorm(0, isigma2)
    for(i in 2:N){
        B[i] <- max(.1, (B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catches[i-1])*exp(u[i]))
        u[i]~dnorm(0, isigma2)
        logcpue[i]~dnorm(log(B[i]*q), itau2)
    }
    B1990 <- B[N]+B[N]*r*(1-B[N]/K)-catches[N]
    P1990 <- B1990/K
}
## temp <- jags(data=data.jags, inits=list(init), param=params.jags, model.file=model.jags,
##             n.chains=1, n.burnin=0, n.iter=5, n.thin=1)
## .jags <- data.frame(temp$BUGSoutput$sims.list)
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
## temp <- stan(fit=model.stan, data=data.stan, iter=5000, chains=1,
##                    warmup=1000, thin=1, init=inits.stan)
rm(dat, data, init,  logk.mean, logk.sigma, logLinf.mean, logLinf.sigma,
   sample.ages, sample.lengths, sample.vbgf, sigma.obs, t0)
## End of Stan
### ------------------------------------------------------------
message("Finished loading ss_logistic models")

temp <- jags(data=data.jags, inits=inits.jags, param=params.jags, model.file=model.jags,
     n.chains=1, n.burnin=1000, n.iter=100*2000+1000, n.thin=100)

.jags <- data.frame(temp$BUGSoutput$sims.list)
mean(.jags$K)
mean(.jags$r)
par(mfrow=c(3,7), mar=c(2,2,2,0))
trash <- lapply(c('K', 'r', 'q', 'B1990', 'P1990', 'sigma2', 'tau2'), function(x)
    hist(.jags[,x], breaks=30, main=x))
trash <- lapply(c('K', 'r', 'q', 'B1990', 'P1990', 'sigma2', 'tau2'), function(x)
    acf(.jags[,x], main=x))
trash <- lapply(c('K', 'r', 'q', 'B1990', 'P1990', 'sigma2', 'tau2'), function(x)
    plot(.jags[,x], type='l',main=x))

par()
plot(r~K, .jags)


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
plot(.jags$r, .jags$K)
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

