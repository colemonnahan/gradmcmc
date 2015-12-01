
message("Generating growth data")
## Global parameters
logLinf.mean <- log(50)
logk.mean <- log(.1)
t0 <- 5
logLinf.sigma <- .1
logk.sigma <- .2
Ntime <- 100
sigma.obs <- .08

## Generate data
sample.vbgf <- function(ages, Linf, k,  t0=5){
    lengths <- Linf*(1-exp(-k*(ages-t0)))
    loglengths <- log(lengths)+ rnorm(n=length(lengths), mean=0, sd=sigma.obs)
    data.frame(ages=ages, loglengths=loglengths)
}
sample.ages <- function(n.ages) {sample(6:Ntime, size=n.ages, replace=FALSE)}
sample.lengths <- function(Nfish, n.ages){
    Linf.vec <- exp(logLinf.mean + rnorm(n=Nfish, 0, sd=logLinf.sigma))
    k.vec <- exp(logk.mean +rnorm(n=Nfish, mean=0, sd=logk.sigma))
    dat <- ldply(1:Nfish, function(i) cbind(fish=i, sample.vbgf(ages=sample.ages(n.ages), Linf=Linf.vec[i], k=k.vec[i])))
    saveRDS(dat, paste0('growth_data_',Nfish,'.RDS'))
    dat
}
## set.seed(5)
## sigma.obs <- .0001
## dat.plot <- sample.lengths(Nfish=50, n.ages=5)
## g <- ggplot(dat.plot, aes(ages, lengths, group=fish)) +
##     geom_point(alpha=.5)
## ggsave('plots/simulated_growth.png', g, width=9, height=5)
## plot(ddply(dat.plot, .(ages), summarize, CV=sd(lengths)/mean(lengths)))
## sigma.obs <- .08
## set.seed(5)
## dat.plot <- sample.lengths(Nfish=50, n.ages=5)
## g <- ggplot(dat.plot, aes(ages, lengths, group=fish))+
##     geom_point(alpha=.5)
## ggsave('plots/observed_growth.png', g, width=9, height=5)

set.seed(5)
dat <- sample.lengths(Nfish=Nfish, n.ages=5)
data <- list(Nfish=Nfish, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
init <- list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf=rep(logLinf.mean, len=Nfish),
                  logk=rep(logk.mean, len=Nfish))

### ------------------------------------------------------------
message("Loading growth models into the workspace")
## JAGS models
data.jags <- data
params.jags <-
    c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma",
      "sigma_obs", "logk", "logLinf")
inits.jags <- list(init)
model.jags <- function(){
    ## hyperpriors
    logLinf_mean~dnorm(log(50),.5)
    logLinf_sigma~dunif(0,.5)
    logk_mean~dnorm(log(.1),.1)
    logk_sigma~dunif(0,.5)
    ## priors
    sigma_obs~dunif(0,.5)
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Nfish){
        logLinf[i]~dnorm(logLinf_mean, pow(logLinf_sigma, -2))
        logk[i]~dnorm(logk_mean, pow(logk_sigma, -2))
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        Linf[i] <- exp(logLinf[fish[i]])
        k[i] <- exp(logk[fish[i]])
        ypred[i] <- Linf[i]*(1-exp(-k[i]*(ages[i]-5)))
        ## ypred[i] <- logLinf[fish[i]]+
        ##     log( (1-exp(-exp(logk[fish[i]])*ages[i])))
        ## Likelihood of data
        loglengths[i]~dnorm(log(ypred[i]), pow(sigma_obs, -2))
    }
}

## temp <- jags(data=data.jags, inits=list(init), param=params.jags, model.file=model.jags,
##      n.chains=1, n.burnin=1000, n.iter=5000, n.thin=1)
## xx <- data.frame(temp$BUGSoutput$sims.list)
## ## acf(xx)
## par(mfrow=c(2,3))
## myre <- function(x,y) (x-y)/x
## myxlim <- c(-.5,.5)
## hist(myre(xx$logLinf_mean,logLinf.mean), xlim=myxlim)
## hist(myre(xx$logLinf_sigma,logLinf.sigma), xlim=myxlim)
## hist(myre(xx$logk_mean,logk.mean), xlim=myxlim)
## hist(myre(xx$logk_sigma,logk.sigma), xlim=myxlim)
## hist(myre(xx$sigma_obs,sigma.obs), xlim=myxlim)
## plot(data$ages, data$loglengths, pch='.')
## a <- 5:Ntime
## lines(a, log(exp(logLinf.mean)*(1-exp(-exp(logk.mean)*(a-5)))))
## lines(a, log(exp(mean(xx$logLinf_mean))*(1-exp(-exp(mean(xx$logk_mean))*(a-5)))))

## End of JAGS
### ------------------------------------------------------------



### ------------------------------------------------------------
## Stan models
data.stan <- data
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
inits.stan <- list(init)
model.stan <- stan(file='growth.stan', data=data.stan, iter=50, chains=1,
                   warmup=10, thin=1, init=inits.stan)
## temp <- stan(fit=model.stan, data=data.stan, iter=5000, chains=1,
##                    warmup=1000, thin=1, init=inits.stan)
## yy <- data.frame(extract(temp, permuted=FALSE)[,1,])
## par(mfrow=c(2,3))
## myre <- function(x,y) (x-y)/x
## myxlim <- c(-.5,.5)
## hist(myre(yy$logLinf_mean,logLinf.mean), xlim=myxlim)
## hist(myre(yy$logLinf_sigma,logLinf.sigma), xlim=myxlim)
## hist(myre(yy$logk_mean,logk.mean), xlim=myxlim)
## hist(myre(yy$logk_sigma,logk.sigma), xlim=myxlim)
## hist(myre(yy$sigma_obs,sigma.obs), xlim=myxlim)
## plot(data$ages, data$loglengths, pch='.')
## a <- 5:Ntime
## lines(a, log(exp(logLinf.mean)*(1-exp(-exp(logk.mean)*(a-5)))))
## lines(a, log(exp(mean(yy$logLinf_mean))*(1-exp(-exp(mean(yy$logk_mean))*(a-5)))))

## acf(xx)
## acf(yy[,1:6])

## par(mfrow=c(2,3))
## qqplot(yy$logLinf_mean, xx$logLinf_mean); abline(0,1)
## qqplot(yy$logLinf_sigma, xx$logLinf_sigma); abline(0,1)
## qqplot(yy$logk_mean, xx$logk_mean); abline(0,1)
## qqplot(yy$logk_sigma, xx$logk_sigma); abline(0,1)
## qqplot(yy$sigma_obs, xx$sigma_obs); abline(0,1)



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
