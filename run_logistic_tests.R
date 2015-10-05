## setwd(main.dir)
setwd(main.dir)
model <- 'logistic'
setwd(model)
compile(paste0(model, '.cpp'))
dyn.load(dynlib(model))
logistic.parameters <- readRDS(results.file(paste0('parameters.', model,'.TMB.RData')))
logistic.data <- readRDS(results.file(paste0('data.', model,'.TMB.RData')))
logr.true <- logistic.parameters$logr
logK.true <- logistic.parameters$logK
logistic.pars <- list(logr=logr.true, logK=logK.true)
logistic.obj <- MakeADFun(data=logistic.data, parameters=logistic.pars, DLL='logistic')
logistic.obj$env$beSilent()
logistic.obj$fn(logistic.pars)
logistic.obj$fn(c(log(.1), log(1000)))
opt <- do.call(optim, logistic.obj)

## explore surface to see how errors are thrown
ngrid <- 35
logr.seq <- seq(log(.05), log(.15), len=ngrid)
logK.seq <- seq(log(4000), log(6000), len=ngrid)
logistic.grid <- ldply(logr.seq, function(r)
    ldply(logK.seq, function(k){
        ## temp.pars <- list(logr=logr.true, logK=logK.true)
        ## temp.pars <- list(logr=r, logK=k)
        ## temp.obj <- MakeADFun(data=logistic.data, parameters=temp.pars, DLL='logistic')
        cbind(r=exp(r), K=exp(k), NLL=logistic.obj$fn(c(r,k)))}))
ggplot(logistic.grid, aes(r,K, size=exp(-NLL), color=exp(-NLL)))+geom_point()

run_hmc <- function(obj, nsim, eps, L, covar=NULL, seed=NULL, ...){
    if(!is.null(seed)) set.seed(seed)
    x1 <- run_mcmc(obj, nsim=nsim, algorithm="HMC", L=L,
                   eps=eps, diagnostic=TRUE, covar=covar)
    minESS <- min(as.vector(effectiveSize(x1$par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=eps, algorithm='hmc', L=L, seed=seed,
                     time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
                     perf=log10(x1$time/minESS))
    return(x2)
}
run_nuts <- function(obj, nsim, covar=NULL, delta, seed=NULL, Madapt, ...){
    if(!is.null(seed)) set.seed(seed)
    x1 <- run_mcmc(obj, nsim=nsim, algorithm="NUTS", diagnostic=TRUE,
                   covar=covar, delta=delta, Madapt=Madapt, ...)
    minESS <- min(as.vector(effectiveSize(x1$par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=delta, algorithm='nuts',
                     seed=seed, time=x1$time, minESS=minESS, acceptance=NA,
                     perf=log10(x1$time/minESS))
    return(x2)
}
run_rwm <- function(obj, nsim, alpha, covar, seed=NULL, ...){
    if(!is.null(seed)) set.seed(seed)
    x1 <- run_mcmc(obj, nsim=nsim, algorithm="RWM", diagnostic=TRUE,
                   covar=covar, alpha=alpha, ...)
    minESS <- min(as.vector(effectiveSize(x1$par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=alpha, algorithm='rwm', seed=seed,
                     time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
                     perf=log10(x1$time/minESS))
    return(x2)
}
run_simulation <- function(model, obj, covar, nsim, seeds.vec, alpha.vec, eps.vec,
                               delta.vec, L.vec){
    Madapt <- min(.1*nsim, 500)
    obj$env$beSilent()
    covar.list <- list(covar, NULL)
    rwm.perf <- ldply(seeds.vec, function(seed)
        ldply(covar.list, function(cov)
            ldply(alpha.vec, function(alpha)
                run_rwm(obj=obj, nsim=nsim, alpha=alpha, seed=seed, covar=cov))))
    hmc.perf <- ldply(seeds.vec, function(seed)
        ldply(covar.list, function(cov)
            ldply(L.vec, function(L)
                ldply(eps.vec, function(eps)
                    run_hmc(obj=obj, nsim=nsim, L=L, eps=eps, seed=seed, covar=cov)))))
    ## nuts.perf <- ldply(seeds.vec, function(seed)
    ##     ldply(covar.list, function(cov)
    ##         ldply(delta.vec, function(delta)
    ##             run_nuts(obj=obj, nsim=nsim, seed=seed, covar=cov,
    ##                      Madapt=Madapt, delta=delta, max_doublings=4))))
    nuts.perf <- NULL
    hmc.perf$algorithm <- paste(hmc.perf$algorithm, hmc.perf$L, sep="_")
    hmc.perf$L <- NULL
    perf <- cbind(model=model, rbind(rwm.perf, hmc.perf, nuts.perf))
    return(perf)
}


logistic.obj$hessian <- TRUE
logistic.opt <- do.call('optim', logistic.obj)
logistic.covar <- solve(logistic.opt$hessian)
nstep <- 10
eps.vec <- seq(.005, .006, len=nstep)
seeds.vec <- c(1,60, 40124,434, 34)[1]
alpha.vec <- seq(.3, 4, len=nstep)
delta.vec <- seq(.1,.4, len=nstep)[1]
L.vec <- c(1,5, 25)
nsim <- 5000
logistic.perf <- run_simulation(model='logistic', obj=logistic.obj, covar=logistic.covar,
               nsim=nsim, seeds.vec=seeds.vec, alpha.vec=alpha.vec,
               eps.vec=eps.vec, delta.vec=delta.vec, L.vec=L.vec)
ggplot(logistic.perf, aes(x=tuning, y=time, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')
ggplot(logistic.perf, aes(x=tuning, y=minESS, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')
ggplot(logistic.perf, aes(x=tuning, y=perf, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')






TMB::find.epsilon(theta=logistic.opt$par, fn=logistic.obj$fn, gr=logistic.obj$gr, eps=.1)
## Run them with and without covariance across gradient of eps
logistic.nuts.list <- logistic.nuts2.list <- logistic.hmc.list <- logistic.hmc2.list <- list()
for(i in 1:length(eps.vec)){
    logistic.hmc.list[[i]] <-
        mcmc(logistic.obj, nsim=nsim, algorithm="HMC", L=10,
             eps=eps.vec[i], params.init=c(logr.true, logK.true),
             diagnostic=TRUE)
    logistic.hmc2.list[[i]] <-
        mcmc(logistic.obj, nsim=nsim, algorithm="HMC", L=10,
             eps=eps.vec[i], params.init=c(logr.true, logK.true),
             diagnostic=TRUE, covar=logistic.covar)
    logistic.nuts.list[[i]] <-
        mcmc(logistic.obj, nsim=nsim, algorithm="NUTS",
             eps=eps.vec[i], params.init=c(logr.true, logK.true),
             diagnostic=TRUE)
    logistic.nuts2.list[[i]] <-
        mcmc(logistic.obj, nsim=nsim, algorithm="NUTS",
             eps=eps.vec[i], params.init=c(logr.true, logK.true),
             diagnostic=TRUE, covar=logistic.covar)
}
## Calculate ESS for each chain, find minimum and then divide by n.calls
## which is approximately the function calls.
temp1 <-laply(logistic.hmc.list, function(x) {log10(x$n.calls/min(as.vector(effectiveSize(x$par))))})
temp2 <-laply(logistic.hmc2.list, function(x) {log10(x$n.calls/min(as.vector(effectiveSize(x$par))))})
logistic.hmc.perf1 <- data.frame(algorithm='hmc', covar=FALSE, eps=eps.vec, perf=temp1)
logistic.hmc.perf2 <- data.frame(algorithm='hmc', covar=TRUE, eps=eps.vec, perf=temp2)
logistic.hmc.perf <- rbind(logistic.hmc.perf1, logistic.hmc.perf2)
ggplot(logistic.hmc.perf, aes(eps, perf, group=covar, color=covar))+ geom_point()
acf(logistic.hmc)
acf(logistic.nuts)
mcmc.nuts2(1, logistic.obj$fn, gr=as.vector(logistic.obj$gr),
          params.init=as.vector(logistic.opt$par), eps=.01)
plot(exp(logistic.mcmc))
points(exp(logr.true), exp(logK.true), pch=16, col='red')
