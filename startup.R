## This assumes the working directory is in the same folder as this file.
## main.dir <- paste0(getwd(), '/')
main.dir <- 'C:/Users/Cole/gradmcmc/'
results.file <- function(file) paste0(main.dir,'results/', file)
plots.file <- function(file) paste0(main.dir,'plots/', file)
figures.file <- function(file) paste0(main.dir,'paper/figures/', file)

#' @param d Dimension of matrix
#' @param covar The covariance for off diagonals. Assumed the same for all
#' pairs of parameters.
#' @return A Hessian matrix (inverse of covariance matrix).
make.covar <- function(d, covar){
    x <- matrix(covar, nrow=d, ncol=d)
    diag(x) <- 1
    x
}
## load libraries and such
library(coda)
library(TMB)
library(ggplot2)
library(plyr)
library(rstan)
library(R2jags)
library(snowfall)
library(reshape2)

## TMB bounding functions, copied from  ADMB
boundpinv <- function(x, min, max){
    -log( (max-min)/(x-min) -1)
}
boundp <- function(x, min, max){
    min + (max-min)/(1+exp(-x))
}

## Get cumulative minimum ESS as a percentage
cum.minESS <- function(df, breaks=5){
    x <- floor(seq(1, dim(df)[1], len=breaks))
    x <- x[x>1]
    y <- lapply(1:length(x), function(i){
       df2 <- df[1:x[i],,, drop=FALSE]
       ess <- data.frame(monitor(df2, warmup=0, print=FALSE))$n_eff
      data.frame(iteration=x[i], pct.ess=100*min(ess)/dim(df2)[1])
 })
    do.call(rbind, y)
}

## wrappers to run chains and return ESS and other metrics
run_hmc <- function(obj, nsim, eps, L, covar=NULL, seed=NULL, diag=FALSE){
    if(!is.null(seed)) set.seed(seed)
    x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="HMC", L=L,
                    eps=eps, diagnostic=TRUE, covar=covar, Madapt=Madapt)
    ## discard the warmup
    par <- x1$par[-(1:Madapt),]
    minESS <- min(as.vector(coda::effectiveSize(par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=eps, algorithm='hmc', L=L, seed=seed,
                     time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
                     perf=log10(x1$time/minESS))
    if(diag) return(x1) else return(x2)
}
run_nuts <- function(obj, nsim, inits=NULL, covar=NULL, delta, seed=NULL,
                     Madapt, diag=FALSE, max_doubling=4){
    if(!is.null(seed)) set.seed(seed)
    x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="NUTS", diagnostic=TRUE, max_doubling=max_doubling,
                   covar=covar, delta=delta, Madapt=Madapt, params.init=inits)
    ## discard the warmup
    par <- x1$par[-(1:Madapt),]
    minESS <- min(as.vector(coda::effectiveSize(par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=delta, algorithm='nuts',
                     seed=seed, time=x1$time, minESS=minESS, acceptance=NA,
                     perf=log10(x1$time/minESS))
    if(diag) return(x1) else return(x2)
}
run_rwm <- function(obj, nsim, alpha, covar, seed=NULL, ...){
    if(!is.null(seed)) set.seed(seed)
    x1 <- TMB::mcmc(obj, nsim=nsim, algorithm="RWM", diagnostic=TRUE,
                   covar=covar, alpha=alpha, ...)
    minESS <- min(as.vector(coda::effectiveSize(x1$par)))
    x2 <- data.frame(covar=!is.null(covar), tuning=alpha, algorithm='rwm', seed=seed,
                     time=x1$time, minESS=minESS, acceptance=mean(x1$accepted),
                     perf=log10(x1$time/minESS))
    return(x2)
}
