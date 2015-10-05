### ------------------------------------------------------------
## This is the central file for running the MCMC efficiency analysis for my
## dissertation.
##
## Started 8/2015

### ------------------------------------------------------------
### Step 0: prepare working space
devtools::install('../adcomp/TMB')
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
covar.temp <- 0
d <- 5
mvn.covar <- make.covar(d=d, covar=covar.temp)
data.TMB <- list(d=d, covar=make.covar(d=d, covar=covar.temp))
parameters.TMB <- list(X=rep(0, times=d))
mvn.obj <- MakeADFun(data=data.TMB, parameters=parameters.TMB, DLL=model)
mvn.obj$env$beSilent()

nsim.seq <- 2^(4:18)
speed.tests <-
    ldply(list(NULL, covar=mvn.covar), function(cov) {
          ldply(nsim.seq, function(x) {
  ## a max_doubling of 5 means 32 steps, so adjusting HMC and RWM
  ## accordingly to make them more comparable for this purpose
  xx <- run_mcmc(obj=mvn.obj, nsim=x, algorithm='NUTS', diag=TRUE, max_doubling=5,
                 delta=.5, Madapt=1, eps=.8, covar=cov)
  x1 <- data.frame(alg='NUTS', nsim=x, time=xx$time, n.calls=xx$n.calls,
                   time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
  xx <- run_mcmc(obj=mvn.obj, nsim=x, algorithm='HMC', diag=TRUE, L=32,
                 eps=.8, covar=cov)
  x2 <- data.frame(alg='HMC', nsim=x, time=xx$time, n.calls=xx$n.calls,
                   time.per.call=xx$time/xx$n.calls, covar=!is.null(cov))
  xx <- run_mcmc(obj=mvn.obj, nsim=32*x, algorithm='RWM', diag=TRUE,
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
xx <- run_mcmc(obj=mvn.obj, nsim=100, algorithm='NUTS', diag=TRUE, max_doubling=5,
               delta=.5, Madapt=1, eps=.8, covar=NULL)
xx <- run_mcmc(obj=mvn.obj, nsim=10, algorithm='RWM', diag=TRUE, alpha=.01,
               covar=NULL)
xx <- run_mcmc(obj=mvn.obj, nsim=10, algorithm='RWM', diag=TRUE, alpha=.01,
               covar=mvn.covar)

Rprof(NULL)
summaryRprof()

x




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


