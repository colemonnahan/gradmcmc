## Generate simulated data
covar <- matrix(cor.temp, nrow=Npar, ncol=Npar)
diag(covar) <- 1
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'
