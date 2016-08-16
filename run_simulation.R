### ------------------------------------------------------------
## This is the central file for running the MCMC efficiency analysis for my
## dissertation.
##
## Started 8/2015

### ------------------------------------------------------------
### Step 0: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'D:/gradmcmc/'
main.dir <- 'C:/Users/Cole/gradmcmc/'
setwd(main.dir)
source("startup.R")
Nout.ind <- 1000
seeds <- c(1:10)
lambda.vec <- NULL
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
## Suppress JAGS and Stan output to file? Useful after debugging to clean
## up console and judge progress.
sink <- TRUE
version$version.string                  # R version 3.2.3
packageVersion('rstan')                 # 2.8.2
packageVersion('R2jags')                # 0.5.7
packageVersion('rjags')                 # 4.4
### End of Step 0.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 1
## Run multivariate normal, empirical and simulated

m <- 'mvnd'
verify <- FALSE
delta <- 0.8
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
## cor is a factor for independent (0) or from wishart (1)
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
verify <- FALSE
delta <- 0.8
Npar <- 5
Nout <- 20000; Nthin <- 1; Nthin.ind <- 10
cor.vec <- c(0, .25, .5, .75, .85, .9, .95)
Npar.vec <- c(2, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 75, 100, 150, 200)
verify <- FALSE
delta <- 0.95
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
verify <- FALSE
delta <- 0.95
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
delta <- 0.95
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
verify <- TRUE
delta <- 0.8
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows; Example 14.5 from Korner-Nievergelt et al
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
verify <- FALSE
delta <- .8
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

### End of Step 1.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 2: Load and prepare result data frames for plotting and tables
source('load_data.R')
source('make_plots.R')
### End of Step 2.
### ------------------------------------------------------------


### ------------------------------------------------------------
## Development code
test <- stan(file='swallows_nc.stan', data=data, init=inits,
             pars=params.jags, iter=2000, chains=1)

sims.stan.nuts <- extract(test, permuted=FALSE)
perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
names.temp <- row.names(perf.stan.nuts)[which(order(perf.stan.nuts$n_eff)<5)]

plot(test, pars=params.jags[1:7])
rstan::traceplot(test, pars=params.jags[1:7])
pairs(test, pars=names.temp)
library(shinystan)
launch_shinystan(test)

test <- stan(file='quantgene_nc.stan', data=data, init=inits,
             pars=pars, iter=2000, chains=1, control=list(adapt_delta=.9))
sims.stan.nuts <- extract(test, permuted=FALSE)
perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
eff <- sort(perf.stan.nuts$n_eff)
barplot(eff[1:50])
names.temp <- row.names(perf.stan.nuts)[which(order(perf.stan.nuts$n_eff)<15)]
plot(test, pars=pars[1:7])
rstan::traceplot(test, pars=names.temp)
pairs(test, pars=c('sa2',names.temp))
acf(sims.stan.nuts[,1, 'sm2'])
library(shinystan)
launch_shinystan(test)

## make sure getting ndivergent right
setwd(main.dir)
setwd('models/mvnd')
Npar <- 5
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'
fit <- stan(file='mvnd.stan', data=data, chains=1, iter=500, warmup=50,
            control=list(adapt_delta=.5, max_treedepth=4))
x <- data.frame(get_sampler_params(fit))
xx <- extract(fit, permuted=FALSE)
sum(x$n_divergent__[-(1:50)])
fit2 <- jags(model.file='mvnd.jags', data=data, inits=inits, para=params.jags,
             n.chains=1, n.iter=5000, n.burnin=2500, n.thin=1)
y <- fit2$BUGSoutput$sims.matrix

test <- run.chains(model='mvnd', seeds=2, Nout=100, delta=.95, data=data,
                   inits=inits, lambda=NULL, params.jags=params.jags,
                   max_treedepth=4)

## explore centering vs noncentering of growth
x <- readRDS('models/growth/sims.ind.RDS')
y <- readRDS('models/growth_nc/sims.ind.RDS')
par(mfrow=c(2,2))
plot(x$logLinf_sigma, x[,'logLinf.1.'])
plot(x$logk_sigma, x[,'logk.1.'])
plot(y$logLinf_sigma, y[,'logLinf_raw.1.'])
plot(y$logk_sigma, y[,'logk_raw.1.'])
par(mfrow=c(2,2))
plot(x$logLinf_sigma, x$logLinf_mean)
plot(x$logk_sigma, x$logk_mean)
plot(y$logLinf_sigma, y$logLinf_mean)
plot(y$logk_sigma, y$logk_mean)

setwd('models/redkite/')
.stan <- readRDS('stan_nuts_diag_e_0.5_1_.RDS')
x <- data.frame(get_sampler_params(.stan))
sum(x$n_divergent__[-(1:5000)])
sum(x$n_divergent__)
.jags <- readRDS('jags_1_.RDS')
sims <- .jags$BUGSoutput$sims.array



## Quick exploration of parallel code
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
Nout.ind=200; Nthin.ind=2000
setwd('models/swallows')
data <- readRDS('data.RDS')
inits <- list(list(a=rep(3.5, len=data$K-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9,
                   fameffphi=rep(0, len=data$nfam),
                   fameffp=rep(0, len=data$nfam),
                   yeareffphi=rep(0, len=4)))
inits <- rep(inits, times=10)
## JAGS complains about environments so easiest is to pass everything as
## absolute numbers.
inits.fn <- function() list(a=rep(3.5, len=18-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9,
                   fameffphi=rep(0, len=72),
                   fameffp=rep(0, len=72),
                   yeareffphi=rep(0, len=4))
params.jags <-
    c('a', 'a1', 'b0', 'b1', 'sigmayearphi', 'sigmaphi', 'sigmap',
      'fameffphi', 'fameffp', 'yeareffphi')

stan.swallows.fit <- stan(file='swallows.stan', par=params.jags, init=inits, data=data,
              iter=20000, thin=10, chains=10, control=list(adapt_delta=.95))
my_sso <- shinystan::launch_shinystan(stan.fit)

xx <- stan(file='swallows.stan', par=params.jags, init=inits, data=data,
              iter=Nout.ind*Nthin.ind, thin=Nthin.ind, chains=10)
my_sso <- shinystan::launch_shinystan(xx)
temp <- extract(xx, permuted=FALSE)
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
perf.stan <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))
sims.stan <- as.data.frame(temp);
names(sims.stan) <- names(xx)

list2env(data, envir=globalenv())
list2env(inits, envir=globalenv())
yy <- jags.parallel(data=data, parameters.to.save=params.jags, inits=inits.fn,
           model.file='swallows.jags', n.chains=10, n.iter=200*2000,
           n.thin=2000)
temp <- yy$BUGSoutput$sims.array
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
sims.jags <- as.data.frame(yy$BUGSoutput$sims.matrix)
perf.jags <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))

perf.platforms <- rbind(cbind(platform='jags',perf.jags),
                        cbind(platform='stan',perf.stan))
perf.platforms <- melt(perf.platforms, c('Rhat', 'n_eff'), id.vars='platform')
plot.model.comparisons(sims.stan, sims.jags, perf.platforms)
