## Took their script and tweaked it to save data and inits
library(R2jags)
inits <- list(list(a=rep(3, len=17), a1=.72, b0=rep(2, len=4),
                   b1=rep(0,len=4),sigmayearphi=.72, sigmaphi=.5,
                   sigmap=.9))
params <- names(inits[[1]])
.jags <- jags(model.file='swallows.jags', data=datax, inits=inits,
              parameters.to.save=params, n.chains=1, n.iter=200)
xx <- as.data.frame(.jags$BUGSoutput$sims.matrix)





##*******************************************************************************
# 14.5  Fledgling Survival of Barn Swallows Using Cormack-Jolly-Seber model in Stan
#*******************************************************************************
#
# Korner-Nievergelt, F.,  T. Roth, S. von Felten,  J. Guelat, B. Almasi and
# P. Korner-Nievergelt. 2015. Bayesian Data Analysis in Ecology using Linear
# Models with R, BUGS and Stan. Elsevier.
#
# Last modification: 2 October 2014
#-------------------------------------------------------------------------------


library(blmeco)
library(rstan)


#-----------------------------------------------------------------------------------
# prepare data
#-----------------------------------------------------------------------------------

data(survival_swallows)

datax <- survival_swallows
str(datax)
inits <- list(list(a=rep(3.5, len=17), a1=.72, b0=rep(0, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9))

saveRDS(datax, 'data.RDS')
saveRDS(inits, 'inits.RDS')
.jags <- jags(data=datax, inits=inits, model.file='swallows.jags',
              n.chains=1, n.iter=1000, n.burnin=500, n.thin=10,
              parameters.to.save=names(inits[[1]]),
              )
xx <- as.data.frame(.jags$BUGSoutput$sims.matrix)
apply(xx, 2, mean)
##-------------------------------------------------------------------------------
# model fit using stan
#-------------------------------------------------------------------------------



# next line needs around 10 min, alternatively you can read in the
# fitted model using the load-function (line below)
mod <- stan(file = "CJS_swallows.stan", data=datax, chains=1,
            iter=1000)
mod <- stan(file = "STAN/CJS_swallows.stan", data=datax, chains=4,
            iter=5000, thin=10, fit=mod)

# load("results/fitted_cjs_swallows.rda")
#save(mod, file="../results/fitted_cjs_swallows.rda")

#-------------------------------------------------------------------------------
# look at results
#-------------------------------------------------------------------------------

           results <- as.data.frame(rstan::extract(mod))
           mod
print(mod, c("a", "a1", "b0", "b1",
             "sigmayearphi", "sigmaphi", "sigmap"), dig=2)


#-------------------------------------------------------------------------------
# look at convergence
#-------------------------------------------------------------------------------

traceplot(mod, "sigmap")
traceplot(mod, "sigmaphi")
traceplot(mod, paste0("a[", 1:8,"]"))

#or alternatively:
library(blmeco)
historyplot(mod, "sigmap")
historyplot(mod, "sigmaphi")
historyplot(mod, "a")
historyplot(mod, "b0")

#-------------------------------------------------------------------------------
# predictive model checking
#-------------------------------------------------------------------------------

# simulate data from the model
modsims <- extract(mod)
nsim <- dim(modsims$phi)[1]  # extract the number of simulations
yrep <- array(dim=c(nsim, datax$I, datax$K))
zrep <- array(dim=c(nsim, datax$I, datax$K))
yrep[,,1] <- 1 # first occasion is 1 for all individuals
zrep[,,1] <- 1
for(j in 1:(datax$K-1)){
  zrep[,,j+1] <- rbinom(nsim*datax$I, size=zrep[,,j], prob=modsims$phi[,,j])
  yrep[,,j+1] <- rbinom(nsim*datax$I, size=zrep[,,j+1], prob=modsims$p[,,j])
}



nobspind <- apply(datax$CH, 1, sum)  # number of observaions per individual
mpfam <- tapply(nobspind, datax$family, mean)
minpfam <- tapply(nobspind, datax$family, min)
maxpfam <- tapply(nobspind, datax$family, max)


repnpind <- apply(yrep, c(1,2), sum)
repmpfam <- matrix(nrow=nsim, ncol=datax$nfam)
repminpfam <- matrix(nrow=nsim, ncol=datax$nfam)
repmaxpfam <- matrix(nrow=nsim, ncol=datax$nfam)

for(f in 1:nsim) {
  repmpfam[f,] <- tapply(repnpind[f,], datax$family, mean)
  repminpfam[f,] <- tapply(repnpind[f,], datax$family, min)
  repmaxpfam[f,] <- tapply(repnpind[f,], datax$family, max)
}

meanrepmpfam <- apply(repmpfam, 2, mean)
meanrepminpfam <- apply(repminpfam, 2, mean)
meanrepmaxpfam <- apply(repmaxpfam, 2, mean)



tiff(filename = "../figures/Figure14-6_predcheck_survival.tif",
     width = 2000, height = 2000, units = "px", pointsize = 45,
     compression = "none")
index <- order(mpfam)
par(mar=c(4,5.5,1,1))
plot(1:72, mpfam[index], xaxt="n", xlab="Family", cex.lab=1.6, cex.axis=1.4,
     type="l", lwd=6, ylab="Number of observations per individual", las=1)
axis(1, at=1:72, labels=c(1:72)[index], cex.axis=1.4)
lines(1:72, minpfam[index], lwd=3)
lines(1:72, maxpfam[index], lwd=3)

lines(1:72, meanrepmpfam[index], col="blue", lwd=6)
lines(1:72, meanrepminpfam[index], col="blue", lwd=3)
lines(1:72, meanrepmaxpfam[index], col="blue", lwd=3)
dev.off()


#-------------------------------------------------------------------------------
# draw inference
#-------------------------------------------------------------------------------


ests <- plogis(apply(modsims$a, 2, mean))
ests.lwr <- plogis(apply(modsims$a, 2, quantile, prob=0.025))
ests.upr <- plogis(apply(modsims$a, 2, quantile, prob=0.975))


# effect of parental care on daily survival
ma <- apply(modsims$a[,1:12], 1, mean) # mean intercept for the first 12 days
newdat <- data.frame(carez=seq(-2, 2, length=100))
b <- c(mean(ma), mean(modsims$a1))
Xmat <- model.matrix(~carez, data=newdat)
newdat$fit <- plogis(Xmat%*%b)
nsim <- nrow(modsims$a)
fitmat <- matrix(ncol=nsim, nrow=nrow(newdat))
for(i in 1:nsim) fitmat[,i] <- plogis(Xmat%*%c(ma[i], modsims$a1[1]))
newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)



# daily survival
tiff(filename = "../figures/Figure14-x_results_survival.tif",
     width = 2000, height = 1000, units = "px", pointsize = 45,
     compression = "none")

par(mfrow=c(1,2), mar=c(5,4.5, 1,1))
plot(1:17, seq(0.5, 1, length=17), type="n", ylab="Daily survival",
     xlab="Day after fleding", las=1)
lines(1:17, ests, lwd=6, col="blue")
lines(1:17, ests.lwr, lwd=3,lty=3, col="blue")
lines(1:17, ests.upr, lwd=3,lty=3, col="blue")



plot(c(-2,2), c(0.6,1), type="n", xlab="Parental care", ylab="Daily survival",
     las=1)
lines(newdat$carez, newdat$fit, lwd=6, col="blue")
lines(newdat$carez, newdat$lwr, lwd=3, lty=3, col="blue")
lines(newdat$carez, newdat$upr, lwd=3, lty=3, col="blue")
dev.off()


# posterior distribution of hypothesis that care has a positive effect:

mean(modsims$a1>0)
nsim/(nsim+1)
