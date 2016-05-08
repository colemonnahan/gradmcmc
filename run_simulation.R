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
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
## cor is a factor for independent (0) or from wishart (1)
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 25, 50, 100, 150, 200)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
verify <- FALSE
delta <- 0.8
Npar <- 5
Nout <- 50000; Nthin <- 1; Nthin.ind <- 10
cor.vec <- c(0, .25, .5, .75, .8, .85, .9, .95, .99, .999)
Npar.vec <- c(2, 50, 100)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 50000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 100)
verify <- FALSE
delta <- 0.8
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic
Nout <- 50000; Nthin <- 1; Nthin.ind <- 500
verify <- FALSE
delta <- 0.8
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
verify <- FALSE
delta <- 0.8
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows; Example 14.5 from Korner-Nievergelt et al
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
verify <- FALSE
delta <- .9
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .9
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

## quantgene; Example 14.5 from Korner-Nievergelt et al
Nout <- 50000; Nthin <- 1; Nthin.ind <- 500
verify <- FALSE
max_treedepth <- 12
m <- 'quantgene_nc'; delta <- .98
source(paste0('models/',m,'/run_model.R'))
m <- 'quantgene'; delta <- .9
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
setwd('models/mvnd')
Npar <- 5
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'
fit <- stan(file='mvnd.stan', data=data, chains=1, iter=50,
            control=list(adapt_delta=.5, max_treedepth=4))
x <- data.frame(get_sampler_params(fit))
sum(x$n_divergent__[-(1:25)])
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
