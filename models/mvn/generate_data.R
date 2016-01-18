## Generate simulated data
df <- max(Npar.vec)
covar <- rWishart(n=1, df=df, Sigma=diag(df))[1:Npar,1:Npar,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'
