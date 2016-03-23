## Generate simulated data
df <- max(Npar.vec)
covar <- rWishart(n=1, df=df, Sigma=diag(df))[1:Npar,1:Npar,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
params.jags <- 'mu'

## cors <- ldply(c(2,5), function(Npar){
##     ldply(c(5,20,200,2000), function(df){
##         ldply(1:5, function(seed){
##             covar <- rWishart(n=1, df=df, Sigma=diag(Npar))[,,1]
##             mycor <- covar/sqrt(diag(covar) %o% diag(covar))
##             diag(mycor) <- 0
##             mycor <- max(abs(mycor))
##             return(cbind(Npar, df, seed, mycor))})})})
## ggplot(cors, aes(df, mycor, group=Npar))+ geom_point()

