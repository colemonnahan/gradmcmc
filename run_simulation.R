### ------------------------------------------------------------
## This is the central file for running the MCMC efficiency analysis for my
## dissertation.
##
## Started 8/2015

### ------------------------------------------------------------
### Step 0: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'C:/Users/Cole/gradmcmc/'
main.dir <- 'D:/gradmcmc/'
setwd(main.dir)
source("startup.R")
Nout.ind <- 1000
seeds <- c(1:6)
lambda.vec <- NULL
delta.vec <- c(.5, .7, .8, .9, .95)
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
### End of Step 0.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 1: Loop through each model and fit empirical and simulated data
m <- 'mvn'
Nout <- 2000; Nthin <- 1; Nthin.ind <- 100
Npar.vec <- c(5, 10, 50, 100, 200, 500, 1000)
source(paste0('models/',m,'/run_model.R'))

m <- 'growth'
Nout <- 2000; Nthin <- 1; Nthin.ind <- 50
Npar.vec <- c(5,10,50, 100, 200, 500, 1000)[-(1:3)]
source(paste0('models/',m,'/run_model.R'))

m <- 'growth2'
Nout <- 2000; Nthin <- 1; Nthin.ind <- 50
Npar.vec <- c(5,10,50, 100, 200, 500, 1000)[1:3]
source(paste0('models/',m,'/run_model.R'))



### End of Step 1.
### ------------------------------------------------------------
### ------------------------------------------------------------
### Step 2: Run simulations. For each model, run across the set of
### algorithms.
### End of Step 2.
### ------------------------------------------------------------



### ------------------------------------------------------------
### Step 3: Load and prepare data
### End of Step 3.
### ------------------------------------------------------------
### ------------------------------------------------------------
### Step 4: Create plots, figures, and tables
### End of Step 2.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
## Play with Neal's funnel from p205 of manual
y <- rnorm(1e3, 0, 3)
x <- rnorm(1e3, 0, sd=exp(y/2))
plot(x, y)
yraw <- rnorm(1e3, 0, 1)
xraw <- rnorm(1e3, 0, 1)
y2 <- 3*yraw
x2 <- exp(y2/2)*xraw
points(x2,y2, col='red')
sigma <- abs(rnorm(1e3, 10, 3))
mu <- 5#rnorm(1e3, 5, 1)
beta <- rnorm(1e3, mu, sigma)
plot(sigma, beta)
beta.raw <- rnorm(1e3,0,1)
points(sigma, mu+sigma*beta.raw, pch=16, col='red')
plot(sigma, beta.raw, pch=16, col='red')

fit.jags <- jags(model.file='mvn.jags', inits=inits.jags, para=params.jags,
                   n.chains=1, data=data.jags)
sims.jags <- fit.jags$BUGSoutput$sims.array
fit.stan <- stan(file='growth2.stan', data=data, iter=5000, chains=1, alg='HMC',
                   thin=1, init=inits, control=list(metric='unit_e'))
sims.stan.hmc <- extract(fit.stan, permuted=FALSE)
minESS.coda <- min(coda::effectiveSize(as.data.frame(sims.stan.hmc[,1,])))

xx <- as.data.frame(get_sampler_params(model.stan))


## source("run_growth.R")
Nyears.vec <- c(20, 50, 100, 200, 500)
L.vec <- c(1,5, 20)
seeds <- c(4,6,1, 48, 10)
source('ss_logistic/run_chains.R')


## Parameters for the short chains, long ones (to verify models are the
## same) are specified in each file and differ.
n.iter1 <- 50000
n.burnin1 <- min(2000,.2*n.iter1)
source("ehook/ehook_run.R")


## Speed tests for Stan algorithms. HMC should scale linearly with int_time
## and NUTS should be similar for average n_leapfrog__ for the warmup and
## sampling period.
seeds <- 1:10
int_time <- seq(1,20, by=2)
Nsim <- 5000
schools_dat <- list(J=8, y=c(28,8,-3,7,-1,1,18,12), sigma=c(15,10,16,11,9,11,10,18))
fit <- stan(file='schools/schools.stan', data=schools_dat, iter=100)
speed.hmc <-
  ldply(seeds, function(seed){
      ldply(int_time, function(L){
          fit2 <- stan(fit=fit, data=schools_dat, iter=Nsim,
                       chains=1, seed=seed, algorithm='HMC', control=list(int_time=L))
          temp <- data.frame(get_sampler_params(fit2))
          data.frame(L=L/temp$stepsize__[Nsim], seed=seed,
                     rstan::get_elapsed_time(fit2))
      })})
speed.nuts <-
    ldply(seeds, function(seed){
        fit2 <- stan(fit=fit, data=schools_dat, iter=Nsim,
                     chains=1, seed=seed, algorithm='NUTS')
        temp <- data.frame(get_sampler_params(fit2))
        data.frame(warmup.nleap=mean(temp[(1:Nsim/2), 'n_leapfrog__']),
                   sample.nleap=mean(temp[-(1:Nsim/2),'n_leapfrog__']),
                   get_elapsed_time(fit2), seed=seed)})
par(mfrow=c(1,2))
plot(warmup~L, data=speed.hmc, ylab="Time", xlab="Steps", main='Warmup',
     ylim=c(0, max(speed.hmc$sample)))
points(warmup~warmup.nleap, data=speed.nuts, col='red', pch=2)
plot(sample~L, data=speed.hmc, ylab="Time", xlab="Steps", main='Sample',
     ylim=c(0, max(speed.hmc$sample)))
points(sample~sample.nleap, data=speed.nuts, col='red', pch=2)
legend('topleft', legend=c("HMC", "NUTS"), col=c('black', 'red'), pch=c(1,2))

## ## Development code
## nsim.seq <- 2^(4:18)
## speed.tests <-
##     ldply(list(NULL, covar=mvn.covar), function(cov) {
##           ldply(nsim.seq, function(x) {
##   ## a max_doubling of 5 means 32 steps, so adjusting HMC and RWM
##   ## accordingly to make them more comparable for this purpose
##   xx <- mcmc(obj=mvn.obj, nsim=x, algorithm='NUTS', diag=TRUE, max_doubling=5,
##                  delta=.5, Madapt=1, eps=.8, covar=cov)
##   x1 <- data.frame(alg='NUTS', nsim=x, time=xx$time, n.calls=xx$n.calls,
##                    time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
##   xx <- mcmc(obj=mvn.obj, nsim=x, algorithm='HMC', diag=TRUE, L=32,
##                  eps=.8, covar=cov)
##   x2 <- data.frame(alg='HMC', nsim=x, time=xx$time, n.calls=xx$n.calls,
##                    time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
##   xx <- mcmc(obj=mvn.obj, nsim=32*x, algorithm='RWM', diag=TRUE,
##                  alpha=.5, covar=cov)
##   x3 <- data.frame(alg='RWM', nsim=x, time=xx$time, n.calls=xx$n.calls,
##                    time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
##   return(rbind(x1,x2,x3))
## })})
## speed.tests$.id <- NULL
## speed.tests$covar <- paste0("covar=",speed.tests$covar)
## saveRDS(speed.tests, results.file('speed.tests.RData'))
## speed.tests <- readRDS(results.file('speed.tests.RData'))
## speed.tests <- reshape2::melt(speed.tests, c('alg', 'nsim', 'covar'))
## ggplot(speed.tests, aes(log2(nsim), log2(value), group=alg, color=alg))+
##     geom_line() + geom_point()+
##     facet_grid(variable~covar, scales='free_y')
## ggsave(plots.file('mvn_speed_comparisons.png'), width=7, height=7)
## speed.tests.normalized <- ddply(speed.tests, .(alg, nsim, variable),
##                                 summarize, ratio=value[2]/value[1])
## ggplot(speed.tests.normalized, aes(log2(nsim), ratio, group=alg, color=alg))+
##     geom_line() + geom_point()+
##     facet_grid(variable~., scales='fixed')
## ## We can give the algorithms information about the covariance of the
## ## posterior to try and improve performance
## covar <- as.matrix(cov(nuts$par))
## rwm2 <- mcmc(obj=obj, nsim=500*8, algorithm='RWM', params.init=opt$par,
##             alpha=1, diagnostic=TRUE, covar=covar)
## rwm2$par <- rwm2$par[seq(1, nrow(rwm2$par), by=8),]
## hmc2 <- mcmc(obj=obj, nsim=500, algorithm='HMC', L=8, params.init=opt$par,
##             diagnostic=TRUE, covar=covar, eps=1)
## nuts2 <- mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
##              diagnostic=TRUE, covar=covar, eps=.1)
## ## See how they compare via ACF
## par(mfrow=c(3,4))
## for(i in 1:4) acf(rwm2$par[,i])
## for(i in 1:4) acf(hmc2$par[,i])
## for(i in 1:4) acf(nuts2$par[,i])
## ## also look at effective size per time
## min(coda::effectiveSize(rwm2$par))/rwm2$time
## min(coda::effectiveSize(hmc2$par))/hmc2$time
## min(coda::effectiveSize(nuts$par))/nuts$time
## ### ------------------------------------------------------------
## ## Try again but without the Laplace approximation to integrate out the
## ## random effects (u).
## obj <- MakeADFun(data=list(x=x, B=B, A=A),
##                  parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
##                  DLL="simple",
##                  silent=TRUE)
## opt <- nlminb(obj$par, obj$fn, obj$gr)
## ## Run RWM and two gradient based algorithms, using adative step size (eps)
## ## for each. Start from the MLE.
## rwm <- mcmc(obj=obj, nsim=5000*8, algorithm='RWM', params.init=opt$par,
##             alpha=.00005, diagnostic=TRUE)
## hmc <- mcmc(obj=obj, nsim=5000, algorithm='HMC', L=8, params.init=opt$par,
##             diagnostic=TRUE)
## nuts <- mcmc(obj=obj, nsim=5000, algorithm='NUTS', params.init=opt$par,
##              diagnostic=TRUE)
## ## Thin it to better approximate the gradient methods
## rwm$par <- rwm$par[seq(1, nrow(rwm$par), by=8),]
## plot(nuts$par[,117])
## plot(rwm$par[,117])
## ## See how they compare via ACF
## par(mfrow=c(3,5))
## for(i in 114:118) acf(rwm$par[,i])
## for(i in 114:118) acf(hmc$par[,i])
## for(i in 114:118) acf(nuts$par[,i])
## barplot(coda::effectiveSize(rwm$par))
## ## also look at effective size per time
## min(coda::effectiveSize(rwm$par))/rwm$time
## min(coda::effectiveSize(hmc$par))/hmc$time
## min(coda::effectiveSize(nuts$par))/nuts$time
## ## We can give the algorithms information about the covariance of the
## ## posterior to try and improve performance
## covar <- as.matrix(cov(nuts$par))
## rwm2 <- mcmc(obj=obj, nsim=500*8, algorithm='RWM', params.init=opt$par,
##             alpha=.1, diagnostic=TRUE, covar=covar)
## rwm2$par <- rwm2$par[seq(1, nrow(rwm2$par), by=8),]
## hmc2 <- mcmc(obj=obj, nsim=500, algorithm='HMC', L=8, params.init=opt$par,
##             diagnostic=TRUE, covar=covar)
## nuts2 <- mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
##              diagnostic=TRUE, covar=covar, eps=0.5)
## ## See how they compare via ACF
## ## See how they compare via ACF
## par(mfrow=c(3,5))
## for(i in 114:118) acf(rwm2$par[,i])
## for(i in 114:118) acf(hmc2$par[,i])
## for(i in 114:118) acf(nuts2$par[,i])
## ## also look at effective size per time
## min(coda::effectiveSize(rwm2$par))/rwm2$time
## min(coda::effectiveSize(hmc2$par))/hmc2$time
## min(coda::effectiveSize(nuts$par))/nuts$time
## x
## ## temporary tests/development code:
## nsim <- 200                              # length of chains
## nburn <- floor(.1*nsim)                 # length of burn in
## algs <- c("HMC", "NUTS", "RWM")
## ## run.simulation(model, alg, covar, nsim=nsim, nburn=nburn)
## xx <- run.model('logistic', 'HMC', nsim=nsim, nburn=nburn, covar=NULL)
## yy <- run.model('logistic', 'NUTS', nsim=nsim, nburn=nburn, covar=NULL)
## res <- rbind(xx$df, yy$df)
## ggplot(res, aes(x=alg, y=perf))+geom_point()+facet_grid(model~.)
## xx <- run.model('mvn', 'HMC', nsim=nsim, nburn=nburn, covar=NULL)
## yy <- run.model('mvn', 'NUTS', nsim=nsim, nburn=nburn, covar=NULL)
## res <- rbind(xx$df, yy$df)
## ggplot(res, aes(x=alg, y=perf))+geom_point()+facet_grid(model~.)
## run.model <- function(model, alg, covar, nsim, nburn, seed=NULL, ...){
##     setwd(model);on.exit(setwd(main.dir))
##     ## Check which software will be used and compile/prep data as needed.
##     if(alg %in% c('HMC', 'NUTS')){
##         ## Use TMB
##         compile(paste0(model, '.cpp'))
##         dyn.load(dynlib(model))
##         data.TMB <- readRDS(results.file(paste0('data.', model,'.TMB.RData')))
##         parameters.TMB <- readRDS(results.file(paste0('parameters.', model,'.TMB.RData')))
##         obj <- MakeADFun(data=data.TMB, parameters=parameters.TMB, DLL=model)
##         obj$hessian <- TRUE
##         obj$env$tracemgc <- FALSE
##         opt <- do.call("optim",obj)
##         opt.covar <- solve(opt$hessian)
##     }
##     if(!is.null(seed)) set.seed(seed)
##     if(alg=='HMC'){
##         time <- system.time(
##             mcmc.out <- TMB::mcmc(obj=obj, nsim=nsim, algorithm='HMC',
##                                   L=15, eps=.015, covar=covar))[['elapsed']]
##         mcmc.out <- mcmc.out[-(1:nburn),]
##     }
##     if(alg=='NUTS'){
##         time <- system.time(
##             mcmc.out <- TMB::mcmc(obj=obj, nsim=nsim,
##                                   algorithm='NUTS',
##                                   eps=.015, covar=covar))[['elapsed']]
##         mcmc.out <- mcmc.out[-(1:nburn),]
##     }
##     mcmc.out <- as.data.frame(mcmc.out)
##     ess <- coda::effectiveSize(mcmc.out)
##     min.ess <- as.numeric(min(ess))
##     perf <- time/min.ess
##     results <- list(df=data.frame(model=model, perf=perf, min.ess=min.ess, time= time, alg=alg,
##                     covar=!is.null(covar)), mcmc.out=mcmc.out)
##     return(results)
## }
## runExample('simple')
## library(microbenchmark)
## obj$env$beSilent()
## microbenchmark(obj$fn(opt$par), obj$gr(opt$par))
