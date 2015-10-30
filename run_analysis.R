### ------------------------------------------------------------
## This is the central file for running the MCMC efficiency analysis for my
## dissertation.
##
## Started 8/2015

### ------------------------------------------------------------
### Step 0: prepare working space
## To reinstall TMB locally, go to that folder and run 'make install' from
## the shell. This is preferred to devtools::install. @@@!! ALSO RESTART R OR
## CHANGES MAY NOT BE THERE !!@@@
source("startup.R")
### End of Step 0.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 1: Run examples to generate data. These are dropped in the resutls
### folder and read in during the simulations.
seeds.list <- list("logistic"=352)
setwd(main.dir)
source('generate.data.R')
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
setwd(main.dir)
model <- 'mvn'
setwd(model)
compile(paste0(model, '.cpp'))
dyn.load(dynlib(model))
d <- 5
set.seed(1)
mvn.covar <- rWishart(n=1, df=d, Sigma=diag(d))[,,1]
data.TMB <- list(d=d, covar=mvn.covar)
parameters.TMB <- list(X=rep(0, times=d))
mvn.obj <- MakeADFun(data=data.TMB, parameters=parameters.TMB, DLL=model)
mvn.obj$env$beSilent()
mvn.opt <- optim(mvn.obj$par, mvn.obj$fn, mvn.obj$gr, hessian=TRUE)
solve(mvn.opt$hessian) - mvn.covar

## Basic tests of algorithms
library(TMB)
mvn.rwm <- mcmc(mvn.obj, nsim=5000, algorithm='RWM', diagnostic=TRUE)
mvn.hmc <- mcmc(mvn.obj, nsim=5000, L=16, algorithm='HMC', eps=NULL, diagnostic=TRUE)
## The trace of eps
with(mvn.hmc, plot(epsbar))
mvn.nuts <- mcmc(mvn.obj, nsim=5000, algorithm='NUTS', eps=NULL, diagnostic=TRUE, max_doublings=4)
with(mvn.nuts, plot(epsbar))
mvn.rwm2 <- mcmc(mvn.obj, nsim=5000, algorithm='RWM', diagnostic=TRUE, covar=mvn.covar)
mvn.hmc2 <- mcmc(mvn.obj, nsim=5000, L=16, algorithm='HMC', eps=NULL, diagnostic=TRUE, covar=mvn.covar)
mvn.nuts2 <- mcmc(mvn.obj, nsim=5000, algorithm='NUTS', eps=NULL, diagnostic=TRUE, max_doublings=4, covar=mvn.covar)


## See how they compare via ACF
par(mfrow=c(3,5))
for(i in 1:5) acf(mvn.rwm$par[,i])
for(i in 1:5) acf(mvn.hmc$par[,i])
for(i in 1:5) acf(mvn.nuts$par[,i])

## See how they compare via ACF
par(mfrow=c(3,5))
for(i in 1:5) acf(mvn.rwm2$par[,i])
for(i in 1:5) acf(mvn.hmc2$par[,i])
for(i in 1:5) acf(mvn.nuts2$par[,i])

par(mfrow=c(3,5))
for(i in 1:5) {qqplot(mvn.rwm$par[,i], mvn.rwm2$par[,i]);abline(a=0,b=1)}
for(i in 1:5) {qqplot(mvn.hmc$par[,i], mvn.hmc2$par[,i]);abline(a=0,b=1)}
for(i in 1:5) {qqplot(mvn.nuts$par[,i], mvn.nuts2$par[,i]);abline(a=0,b=1)}


nsim.seq <- 2^(4:18)
speed.tests <-
    ldply(list(NULL, covar=mvn.covar), function(cov) {
          ldply(nsim.seq, function(x) {
  ## a max_doubling of 5 means 32 steps, so adjusting HMC and RWM
  ## accordingly to make them more comparable for this purpose
  xx <- mcmc(obj=mvn.obj, nsim=x, algorithm='NUTS', diag=TRUE, max_doubling=5,
                 delta=.5, Madapt=1, eps=.8, covar=cov)
  x1 <- data.frame(alg='NUTS', nsim=x, time=xx$time, n.calls=xx$n.calls,
                   time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
  xx <- mcmc(obj=mvn.obj, nsim=x, algorithm='HMC', diag=TRUE, L=32,
                 eps=.8, covar=cov)
  x2 <- data.frame(alg='HMC', nsim=x, time=xx$time, n.calls=xx$n.calls,
                   time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
  xx <- mcmc(obj=mvn.obj, nsim=32*x, algorithm='RWM', diag=TRUE,
                 alpha=.5, covar=cov)
  x3 <- data.frame(alg='RWM', nsim=x, time=xx$time, n.calls=xx$n.calls,
                   time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
  return(rbind(x1,x2,x3))
})})
speed.tests$.id <- NULL
speed.tests$covar <- paste0("covar=",speed.tests$covar)
saveRDS(speed.tests, results.file('speed.tests.RData'))
speed.tests <- readRDS(results.file('speed.tests.RData'))
speed.tests <- reshape2::melt(speed.tests, c('alg', 'nsim', 'covar'))
ggplot(speed.tests, aes(log2(nsim), log2(value), group=alg, color=alg))+
    geom_line() + geom_point()+
    facet_grid(variable~covar, scales='free_y')
ggsave(plots.file('mvn_speed_comparisons.png'), width=7, height=7)

speed.tests.normalized <- ddply(speed.tests, .(alg, nsim, variable),
                                summarize, ratio=value[2]/value[1])
ggplot(speed.tests.normalized, aes(log2(nsim), ratio, group=alg, color=alg))+
    geom_line() + geom_point()+
    facet_grid(variable~., scales='fixed')


## Profiling over NUTS since it is so slow
Rprof()
xx <- mcmc(obj=mvn.obj, nsim=100, algorithm='NUTS', diag=TRUE, max_doubling=5,
               delta=.5, Madapt=1, eps=.8, covar=NULL)
xx <- mcmc(obj=mvn.obj, nsim=10, algorithm='RWM', diag=TRUE, alpha=.01,
               covar=NULL)
xx <- mcmc(obj=mvn.obj, nsim=10, algorithm='RWM', diag=TRUE, alpha=.01,
               covar=mvn.covar)

Rprof(NULL)
summaryRprof()

x


## We can give the algorithms information about the covariance of the
## posterior to try and improve performance
covar <- as.matrix(cov(nuts$par))
rwm2 <- mcmc(obj=obj, nsim=500*8, algorithm='RWM', params.init=opt$par,
            alpha=1, diagnostic=TRUE, covar=covar)
rwm2$par <- rwm2$par[seq(1, nrow(rwm2$par), by=8),]
hmc2 <- mcmc(obj=obj, nsim=500, algorithm='HMC', L=8, params.init=opt$par,
            diagnostic=TRUE, covar=covar, eps=1)
nuts2 <- mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
             diagnostic=TRUE, covar=covar, eps=.1)

## See how they compare via ACF
par(mfrow=c(3,4))
for(i in 1:4) acf(rwm2$par[,i])
for(i in 1:4) acf(hmc2$par[,i])
for(i in 1:4) acf(nuts2$par[,i])

## also look at effective size per time
min(coda::effectiveSize(rwm2$par))/rwm2$time
min(coda::effectiveSize(hmc2$par))/hmc2$time
min(coda::effectiveSize(nuts$par))/nuts$time


### ------------------------------------------------------------
## Try again but without the Laplace approximation to integrate out the
## random effects (u).
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
                 DLL="simple",
                 silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Run RWM and two gradient based algorithms, using adative step size (eps)
## for each. Start from the MLE.
rwm <- mcmc(obj=obj, nsim=5000*8, algorithm='RWM', params.init=opt$par,
            alpha=.00005, diagnostic=TRUE)
hmc <- mcmc(obj=obj, nsim=5000, algorithm='HMC', L=8, params.init=opt$par,
            diagnostic=TRUE)
nuts <- mcmc(obj=obj, nsim=5000, algorithm='NUTS', params.init=opt$par,
             diagnostic=TRUE)
## Thin it to better approximate the gradient methods
rwm$par <- rwm$par[seq(1, nrow(rwm$par), by=8),]
plot(nuts$par[,117])
plot(rwm$par[,117])


## See how they compare via ACF
par(mfrow=c(3,5))
for(i in 114:118) acf(rwm$par[,i])
for(i in 114:118) acf(hmc$par[,i])
for(i in 114:118) acf(nuts$par[,i])

barplot(coda::effectiveSize(rwm$par))

## also look at effective size per time
min(coda::effectiveSize(rwm$par))/rwm$time
min(coda::effectiveSize(hmc$par))/hmc$time
min(coda::effectiveSize(nuts$par))/nuts$time

## We can give the algorithms information about the covariance of the
## posterior to try and improve performance
covar <- as.matrix(cov(nuts$par))
rwm2 <- mcmc(obj=obj, nsim=500*8, algorithm='RWM', params.init=opt$par,
            alpha=.1, diagnostic=TRUE, covar=covar)
rwm2$par <- rwm2$par[seq(1, nrow(rwm2$par), by=8),]
hmc2 <- mcmc(obj=obj, nsim=500, algorithm='HMC', L=8, params.init=opt$par,
            diagnostic=TRUE, covar=covar)
nuts2 <- mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
             diagnostic=TRUE, covar=covar, eps=0.5)

## See how they compare via ACF
## See how they compare via ACF
par(mfrow=c(3,5))
for(i in 114:118) acf(rwm2$par[,i])
for(i in 114:118) acf(hmc2$par[,i])
for(i in 114:118) acf(nuts2$par[,i])

## also look at effective size per time
min(coda::effectiveSize(rwm2$par))/rwm2$time
min(coda::effectiveSize(hmc2$par))/hmc2$time
min(coda::effectiveSize(nuts$par))/nuts$time




x
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




