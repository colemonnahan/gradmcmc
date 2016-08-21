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
Npar.vec <- c(1:10, floor(10^(seq(1.25, 3, by=.25))))
verify <- FALSE
delta <- 0.8
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
delta <- 0.8
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

## Wildflower; from Bolker et al 2010
Nout <- 1000; Nthin <- 1; Nthin.ind <- 500
verify <- FALSE
delta <- .8
m <- 'wildflower_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'wildflower'
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


## Special runs of growth model to show divergences and effect of
## noncentering
logLinf.mean <- log(50); logk.mean <- log(.1); t0 <- 6; logLinf.sigma <- .1; logk.sigma <- .2; Ntime <- 40; sigma.obs <- .1
Nfish <- 10
set.seed(115)
dat <- sample.lengths(Nfish=Nfish, n.ages=5, logLinf.mean=logLinf.mean,
                           logLinf.sigma=logLinf.sigma, logk.mean=logk.mean,
                           logk.sigma=logk.sigma, sigma.obs=sigma.obs, t0=t0)
g <- ggplot(dat, aes(ages, loglengths, group=fish, color=fish)) +
    geom_point(alpha=.5) + geom_line()
data <- list(Nfish=Nfish, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
inits <- list(list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf=rep(logLinf.mean, len=Nfish),
                  logk=rep(logk.mean, len=Nfish), delta=1))
pars <-
  c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma", "logk", "logLinf",
    "sigma_obs", "delta")
iter <- 20000
temp <- stan(file='models/growth/growth.stan', data=data, pars=pars,
             chains=1, iter=1)
growth.fit <- stan(fit=temp, data=data, pars=pars, init=inits, iter=iter,
                   chains=1, control=list(adapt_delta=.9))
growth.sims <- data.frame(extract(growth.fit, permuted=FALSE)[,1,])
growth.divergent <- data.frame(get_sampler_params(growth.fit, inc_warmup=FALSE))
## Noncentered model
inits <- list(list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf_raw=rep(0, len=Nfish),
                  logk_raw=rep(0, len=Nfish), delta=1))
pars <-
    c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma",
      "sigma_obs", "logk_raw","logLinf_raw", "delta")
tempnc <- stan(file='models/growth_nc/growth_nc.stan', data=data, pars=pars,
             chains=1, iter=1)
growth_nc.fit <- stan(fit=tempnc, data=data, pars=pars, init=inits,
             iter=iter, chains=1, control=list(adapt_delta=.9))
growth_nc.sims <- data.frame(extract(growth_nc.fit, permuted=FALSE)[,1,])
growth_nc.divergent <- data.frame(get_sampler_params(growth_nc.fit, inc_warmup=FALSE))
xx <- data.frame(sigma=growth.sims$logLinf_sigma, tau=growth.sims$logLinf.1.-growth.sims$logLinf_mean,
            divergent=growth.divergent$divergent__)
yy <- data.frame(sigma=growth_nc.sims$logLinf_sigma, tau=growth_nc.sims$logLinf_raw.1.,
            divergent=growth_nc.divergent$divergent__)
## Save to results and then plot them elsewhere
saveRDS(xx, file='results/growth_divergences.RDS')
saveRDS(yy, file='results/growth_nc_divergences.RDS')
### ENd of growth divergence
