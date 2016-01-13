## Fit the empirical data

## Data
Npar <- 50
Niter <- 100000
Nwarmup <- Niter/2
Nthin <- 200
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))

## Inits
inits <- list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2)
inits.jags <- inits.stan <- list(inits)
data.jags <- data.stan <- data
### ------------------------------------------------------------
### First get independent samples to verify the models
## JAGS. Will be fit once since no tuning parameters
params.jags <- c('mu')
model.jags <- 'mvn.jags'
fit.jags <- jags(data=data.jags, inits=inits.jags, param=params.jags,
             model.file=model.jags, n.chains=1, n.burnin=Nwarmup, n.iter=Niter,
             n.thin=Nthin)
sims.jags <- data.frame(fit.jags$BUGSoutput$sims.matrix)
array.jags <- fit.jags$BUGSoutput$sims.array

## Stan. Will be fit across tuning parameters delta and mass matrix
model.stan <- 'mvn.stan'
fit.stan <- stan(file=model.stan, data=data.stan, iter=Niter, chains=1,
                   warmup=Nwarmup, thin=Nthin, init=inits.stan)
sims.stan <- data.frame(extract(fit.stan, permuted=FALSE)[,1,])

## Verify -- LEFT OFF HERE
jags.ess <- data.frame(rstan::monitor(sims=fit.jags$BUGSoutput$sims.array, warmup=0, print=FALSE, probs=.5))$n_eff
array.stan <- array(as.matrix(sims.stan), dim=c(nrow(sims.stan),1, ncol(sims.stan)))
stan.ess <- data.frame(rstan::monitor(sims=array.stan, warmup=0, print=FALSE, probs=.5))$n_eff
rbind(jags.ess, stan.ess)
plot.model.comparisons(sims.stan, sims.jags)
## End of independent draws
### ------------------------------------------------------------


## Verify the models are the same (coding errors)

plot.model.comparisons <- function(sims.stan, sims.jags){
    names(sims.stan) <- gsub('\\.', '', x=names(sims.stan))
    names(sims.jags) <- gsub('\\.', '', x=names(sims.jags))
    sims.stan$lp__ <- sims.jags$deviance <- NULL
    par.names <- names(sims.stan)
    sims.jags <- sims.jags[,par.names]
    qq <- ldply(par.names, function(i){
                    temp <- as.data.frame(qqplot(sims.jags[,i], sims.stan[,i], plot.it=FALSE))
                    return(cbind(i,temp))
                })
    g <- ggplot(qq, aes(x,y))+geom_point(alpha=.5) + geom_abline(slope=1, col='red') + facet_wrap('i')
    ggsave('plots/model_comparison.png', width=9, height=5)

}
