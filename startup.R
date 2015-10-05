## This assumes the working directory is in the same folder as this file.
main.dir <- paste0(getwd(), '/')
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
library(TMB)
library(ggplot2)
library(plyr)
library(coda)
## library(snowfall)
## cores <- 4
## sfInit(parallel=TRUE, cpus=cores, type='SOCK')
## sfLibrary(TMB)
## wrapper <- function(eps){
##     run_mcmc(mvn.obj, nsim=nsim, algorithm="HMC", L=L,
##              eps=eps, params.init=mvn.opt$par,
##              diagnostic=TRUE)
## }
## sfExport('run_mcmc', namespace='TMB')
## mvn.hmc.list <- sfLapply(eps.vec, function(x)  run_mcmc(mvn.obj, nsim=nsim, algorithm="HMC", L=L,
##              eps=x, params.init=mvn.opt$par,
##              diagnostic=TRUE))
