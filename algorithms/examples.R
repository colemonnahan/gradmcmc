### This file quickly demonstrates how to use the static HMC and NUTS
### algorithms provided in the mcmc.R file in this folder. Here we use
### simple multivariate normal model with analytical gradients, avoiding
### the need for automatic differentiation. See function delcarations for
### argument definitions.

library(plyr)
library(mvtnorm)
source('mcmc.R')
col.label <- gray(.3)
col.border <- gray(.5)
col.tick <- gray(.4)
col.contour <- gray(.5)
lty.contour <- 2
lwd.contour <- .75


## 2d parabola: y=x1^2+x2^2 with static HMC leapfrog trajectories, with and
## without mass matrix
covar <- matrix(c(1,-.95,-.95,1), 2)
fn <- function(x) -t(x)%*%solve(covar)%*%x
gr <- function(x) as.vector(-solve(covar)%*%x)
ngrid <- 100
x1.seq <- x2.seq <- seq(-3, 3, len=ngrid)
mvn.grid <- ldply(x1.seq, function(x1)
    ldply(x2.seq, function(x2){
        cbind(x1=x1, x2=x2, NLL=-dmvnorm(x=c(x1,x2), sigma=covar, log=TRUE))}))
mvn.grid$NLL[is.nan(mvn.grid$NLL)] <- NA
z.mvn <- matrix(exp(-mvn.grid$NLL), nrow=ngrid, ncol=ngrid, byrow=TRUE)

L.big <- 50
eps.big <- .1
## without mass matrix
mvn.leap1 <- leapfrog(L=L.big, fn=fn, gr=gr, eps=eps.big,
                      init=c(2.3, -2), r.cur=c(-1,.1))
## with it
mvn.leap2 <- leapfrog(L=L.big, fn=fn, gr=gr, eps=eps.big, covar=covar,
                      init=c(2.3, -2), r.cur=c(-1,.1))
mvn.leap.list <- list(mvn.leap1, mvn.leap2)
nlevels <- 8
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
a <- mvn.leap1$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05)
points(a[1,1], a[1,2], col='red', pch=16)
a <- mvn.leap2$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05, col='blue')
points(a[1,1], a[1,2], col='red', pch=16)
mtext(expression(x[1]), 1, line=2.5)
mtext(expression(x[2]), 2, line=2.5)


## Run NUTS with and without mass matrix
nuts.nocovar <- nuts(nsim=500, fn=fn, gr=gr, params.init=c(2.3, -2), delta=.8,
           diagnostic=TRUE)
nuts.covar <- nuts(nsim=500, fn=fn, gr=gr, params.init=c(2.3, -2), delta=.8,
           diagnostic=TRUE, covar=covar)
## The adaptation of the step size
plot(nuts.nocovar$epsbar, ylab='Step size', type='l')
lines(nuts.covar$epsbar, col='blue')
par(mfrow=c(1,2))
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
points(nuts.nocovar$par)
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
points(nuts.covar$par)
