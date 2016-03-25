## Generate simulated data
covar <- matrix(cor.temp, nrow=Npar, ncol=Npar)
diag(covar) <- 1
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- lapply(1:length(seeds), function(i)
  list(mu=mvrnorm(n=1, mu=rep(0, len=Npar), Sigma=covar)))
params.jags <- 'mu'
