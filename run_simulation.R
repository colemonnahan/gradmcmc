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
delta.vec <- .8 #c(.5, .7, .8, .9, .95)
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
## Suppress JAGS and Stan output? Useful after debugging to clean up
## console and judge progress.
sink <- TRUE
.call('get_version', package='rjags')   # JAGS version 4.2.0
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
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
## cor is a factor for independent (0) or from wishart (1)
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 25, 50, 100, 150, 200)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
verify <- FALSE
Npar <- 5
Nout <- 50000; Nthin <- 1; Nthin.ind <- 10
cor.vec <- c(0, .25, .5, .75, .8, .85, .9, .95, .99, .999)
Npar.vec <- c(2, 50, 100)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 100000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 100)
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic
Nout <- 100000; Nthin <- 1; Nthin.ind <- 500
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows; Example 14.5 from Korner-Nievergelt et al
Nout <- 200000; Nthin <- 1; Nthin.ind <- 100
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

## quantgene; Example 14.5 from Korner-Nievergelt et al
Nout <- 2000; Nthin <- 1; Nthin.ind <- 500
m <- 'quantgene_nc'
source(paste0('models/',m,'/run_model.R'))
## m <- 'quantgene'
## source(paste0('models/',m,'/run_model.R'))
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


test <- run.chains(model=m, seeds=1, Nout=100, delta=.9, data=data,
                   inits=inits, lambda=NULL, params.jags=params.jags)
fit <- stan(file='mvnd.stan', data=data, chains=1, iter=50,
            control=list(adapt_delta=.9, max_treedepth=4))
x <- data.frame(get_sampler_params(fit))
