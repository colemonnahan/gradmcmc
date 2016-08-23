### This file quickly demonstrates how to use the static HMC and NUTS
### algorithms provided in the mcmc.R file in this folder. Here we use
### simple multivariate normal model with analytical gradients, avoiding
### the need for automatic differentiation. See function delcarations for
### argument definitions. These functions do not adapt the mass matrix,
### rather the user specifies it.

library(plyr)
library(mvtnorm)
## This file contains a function to make a leapfrog trajectory
source('mcmc.R')

## 2d multivariate normal with static HMC leapfrog trajectories, with and
## without mass matrix
covar <- matrix(c(1,-.95,-.95,1), 2)
covar.inv <- solve(covar)
## Analytical log-densities and gradients, not optimized so really slow
fn <- function(x) -as.vector(dmvnorm(as.vector(x), sigma=covar, log=TRUE))
gr <- function(x) as.vector(covar.inv%*%x)
## Create a fake minimal TMB object
obj <- list(fn=fn, gr=gr)
obj$env$beSilent <- function() NULL
obj$par <- c(1,1)

## Make a grid of points for a contour
ngrid <- 100
x1.seq <- x2.seq <- seq(-3, 3, len=ngrid)
mvn.grid <- ldply(x1.seq, function(x1)
    ldply(x2.seq, function(x2){
        cbind(x1=x1, x2=x2, NLL=-fn(x=c(x1,x2)))}))
## mvn.grid$NLL[is.nan(mvn.grid$NLL)] <- NA
z.mvn <- matrix(exp(mvn.grid$NLL), nrow=ngrid, ncol=ngrid, byrow=TRUE)

## without mass matrix (covar)
traj1 <-
  leapfrog(L=25, fn=fn, gr=gr, eps=.2, init=c(2.3, -2), r.cur=c(-2,.1))
## with it
traj2 <-
  leapfrog(L=5, fn=fn, gr=gr, eps=1, covar=covar, init=c(2.3, -2), r.cur=c(-2,.1))

## Run NUTS with and without mass matrix (covar), using the TMB
## functions. Then plot some key quantities
nsim <- 200
nuts.nocovar <-
  run_mcmc(obj=obj, nsim=nsim, algorithm='NUTS', delta=0.8,
           max_doublings=8, diag=TRUE)
nuts.covar <-
  run_mcmc(obj=obj, nsim=nsim, algorithm='NUTS', delta=0.8,
           max_doublings=8, diag=TRUE, covar=covar)

## Make exploratory plots of trajectories
par(mfrow=c(2,2))
nlevels <- 10
time2 <- 1*(0:5)
time1 <- .2*(0:25)
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=gray(.5), lty=3,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), xlab='x1', ylab='x2')
a <- traj1$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05)
points(a[1,1], a[1,2], col='red', pch=16)
a <- traj2$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05, col='blue')
points(a[1,1], a[1,2], col='red', pch=16)
## K=kintetic energy, P=potential
with(traj1, plot(time1, K, type='b', ylim=range(c(K,P)), ylab='Energy', xlab='Time'))
lines(time1, traj1$P, type='b', col='blue')
with(traj2, plot(time2,K, type='b', ylim=range(c(K,P)), ylab='Energy', xlab='Time'))
lines(time2, traj2$P, type='b', col='blue')
with(traj1, plot(time1[-1],accept, type='l', ylim=c(.5,1.01), ylab='Acceptance Ratio', xlab='Time'))
with(traj2, lines(time2[-1],accept, col='blue'))

## Plots of NUTS
par(mfrow=c(2,2))
plot(nuts.nocovar$epsbar, ylab='Step size', xlab='Warmup Iterations', type='l')
lines(nuts.covar$epsbar, col='blue')
plot(.8-nuts.nocovar$Hbar, ylab='Acceptance Rate', xlab='Warmup Iterations', type='l', ylim=c(.5,1))
lines(.8-nuts.covar$Hbar, col='blue')
abline(h=.8)
plot((nsim/2+1):nsim, log2(nuts.nocovar$steps.taken[-(1:(nsim/2))]), ylab='Tree Depth',
     xlab='Sampling Iterations', type='l', ylim=c(0,8))
lines((nsim/2+1):nsim, log2(nuts.covar$steps.taken[-(1:(nsim/2))]), col='blue')
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=gray(.5), lty=3,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), xlab='x1', ylab='x2')
points(nuts.nocovar$par)
points(nuts.covar$par, col='blue')
