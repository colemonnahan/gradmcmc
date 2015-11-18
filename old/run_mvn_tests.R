## Run experiments on MVN with different correlations and dimensions
setwd(main.dir)
model <- 'mvn'
setwd(model)
compile(paste0(model, '.cpp'))
dyn.load(dynlib(model))
covar.temp <- .8
d <- 2
mvn.covar <- make.covar(d=d, covar=covar.temp)
data.TMB <- list(d=d, covar=make.covar(d=d, covar=covar.temp))
parameters.TMB <- list(X=rep(0, times=d))
mvn.obj <- MakeADFun(data=data.TMB, parameters=parameters.TMB, DLL=model)
mvn.obj$env$beSilent()

## mvn.fn <- function(x) -mvn.obj$fn(x)
## mvn.gr <- function(x) -as.vector(mvn.obj$gr(x))
## mvn.samples <- mvtnorm::rmvnorm(n=1000, sigma=mvn.covar)
## plot(mvn.samples)
## rm(theta.trajectory)
## .buildtree(theta=c(3,3), theta0=c(3,3),  r=c(0,1), r0=c(0,1),  u=100, v=-1,
##            j=8, eps=.01, fn=mvn.fn, gr=mvn.gr)
## z <- nrow(theta.trajectory); x <- theta.trajectory[,1]; y <- theta.trajectory[,2]
## arrows(x[-z], y[-z], x[-1], y[-1], col='red', pch=16, cex=.5, length=.05)

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
run_mvn_simulation <- function(d, covar.temp, nsim, seeds.vec, alpha.vec, eps.vec,
                               delta.vec, L.vec){
    print(d); print(covar.temp)
    Madapt <- min(.1*nsim, 500)
    mvn.covar <- make.covar(d=d, covar=covar.temp)
    data.temp <- list(d=d, covar=mvn.covar)
    mvn.obj <- MakeADFun(data=data.temp, parameters=list(X=rep(0, times=d)), DLL='mvn')
    mvn.obj$env$beSilent()
    covar.list <- list(mvn.covar, NULL)
    mvn.rwm.perf <- ldply(seeds.vec, function(seed)
        ldply(covar.list, function(cov)
            ldply(alpha.vec, function(x)
                run_rwm(obj=mvn.obj, nsim=nsim, alpha=x, seed=seed, covar=cov))))
    mvn.hmc.perf <- ldply(seeds.vec, function(seed)
        ldply(covar.list, function(cov)
            ldply(L.vec, function(L)
                ldply(eps.vec, function(x)
                    run_hmc(obj=mvn.obj, nsim=nsim, L=L, eps=x, seed=seed, covar=cov)))))
    mvn.nuts.perf <- ldply(seeds.vec, function(seed)
        ldply(covar.list, function(cov)
            ldply(delta.vec, function(x)
                run_nuts(obj=mvn.obj, nsim=nsim, seed=seed, covar=cov,
                         Madapt=Madapt, delta=x, max_doublings=4))))
    mvn.hmc.perf$algorithm <- paste(mvn.hmc.perf$algorithm, mvn.hmc.perf$L, sep="_")
    mvn.hmc.perf$L <- NULL
    perf <- cbind(model=paste('mvn', d, covar.temp, sep='_'),
                  rbind(mvn.rwm.perf, mvn.hmc.perf, mvn.nuts.perf))
    return(perf)
}


delta.vec <- seq(.1,.9, len=15)
seeds.vec <- 1
nsim <- 50000
Madapt <- .1*nsim
covar.list <- list(mvn.covar, NULL)
mvn.nuts.perf <-
    ldply(seeds.vec, function(seed)
          ldply(covar.list, function(cov)
                ldply(delta.vec, function(x)
                      run_nuts(obj=mvn.obj, nsim=nsim, seed=seed, covar=cov,
                               Madapt=Madapt, delta=x, max_doublings=6))))
saveRDS(mvn.nuts.perf, results.file('mvn.nuts.perf.RDS'))
mvn.nuts.perf$acceptance <- NULL
mvn.nuts.perf.long <-
    reshape2::melt(mvn.nuts.perf, c('covar', 'tuning', 'algorithm', 'seed'))
ggplot(mvn.nuts.perf.long, aes(tuning, value, group=covar, color=covar)) +
    geom_line() + facet_wrap('variable', nrow=3, scales='free')
ggsave(plots.file('mvn_nuts_perf.png'), width=8, height=6)

## Run preliminary tests on some algorithms
nsim <- 100000
nseq <- 15
seeds.vec <- c(1,60, 40124,434, 34)[1]
## alpha.vec <- seq(.3, 4, len=nseq)
delta.vec <- seq(.01,.99, len=nseq)
alpha.vec <- exp(seq(log(.001), log(10), len=nseq))
eps.vec <- seq(.05, 2, len=nseq)
L.vec <- c(1,5, 15)
perf <- ldply(c(2,20,100), function(d)
    ldply(c(0, .5, .98), function(cov)
        run_mvn_simulation(
            d=d, covar.temp=cov, nsim=nsim,
            seeds.vec=seeds.vec, eps.vec=eps.vec,
            alpha.vec=alpha.vec, delta.vec=delta.vec, L.vec=L.vec)))
saveRDS(perf, file=results.file('mvn.perf.RData'))
mvn.perf <- readRDS(results.file('mvn.perf.RData'))
ggplot(mvn.perf, aes(x=tuning, y=time, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')
ggplot(mvn.perf, aes(x=tuning, y=minESS, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')
ggplot(mvn.perf, aes(x=tuning, y=perf, group=covar, color=covar))+geom_line()+ geom_point(size=.5)+
    facet_grid(model~algorithm, scales='free')
ggsave(plots.file('mvn_perf.png'), width=12, height=9)


## ggplot(mvn.nuts.perf, aes(tuning, minESS, group=covar, color=covar))+ geom_point()
## ggplot(mvn.nuts.perf, aes(tuning, perf, group=covar, color=covar))+ geom_point(apha=.5)
## ggplot(mvn.rwm.perf, aes(tuning, time, group=covar, color=covar))+ geom_point()
## ggplot(mvn.rwm.perf, aes(tuning, minESS, group=covar, color=covar))+ geom_point()
## ggplot(mvn.rwm.perf, aes(tuning, acceptance, group=covar, color=covar))+
##     geom_point() + ylim(0,1)
## ggplot(mvn.rwm.perf, aes(tuning, perf, group=covar, color=covar))+ geom_point(apha=.5)
## ggplot(mvn.hmc.perf, aes(tuning, time, group=covar, color=covar))+ geom_point() +
##     facet_grid(.~L)
## ggplot(mvn.hmc.perf, aes(tuning, minESS, group=covar, color=covar))+ geom_point() +
##     facet_grid(.~L)
## ggplot(mvn.hmc.perf, aes(tuning, acceptance, group=covar, color=covar))+ geom_point() +
##     facet_grid(.~L)
## ggplot(mvn.hmc.perf, aes(tuning, perf, group=covar, color=covar))+ geom_point(apha=.5) +
##     facet_grid(.~L)
## ## extract out the best chain for comparing between algorithms

setwd(main.dir)
