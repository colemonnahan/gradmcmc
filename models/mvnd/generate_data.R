## Generate simulated data
df <- max(Npar.vec)
if(j==0){
    covar <- diag(Npar)
} else if(j==1){
    covar <- rWishart(n=1, df=df, Sigma=diag(df))[1:Npar,1:Npar,1]
} else {
    stop("Invalid cor value in mvnd")
}
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- lapply(1:length(seeds), function(i)
  list(mu=mvrnorm(n=1, mu=rep(0, len=Npar), Sigma=covar)))
params.jags <- 'mu'
### old tests for how to get correlations along a gradient, for now using a
### two step approach (above)
## cors <- ldply(c(2,5), function(Npar){
##     ldply(c(5,20,200,2000), function(df){
##         ldply(1:5, function(seed){
##             covar <- rWishart(n=1, df=df, Sigma=diag(Npar))[,,1]
##             mycor <- covar/sqrt(diag(covar) %o% diag(covar))
##             diag(mycor) <- 0
##             mycor <- max(abs(mycor))
##             return(cbind(Npar, df, seed, mycor))})})})
## ggplot(cors, aes(df, mycor, group=Npar))+ geom_point()

