## A script to prepare data and do some quick checks that the models are
## working correctly.

library(blmeco)
library(rstan)
library(R2jags)
library(plyr)
#-----------------------------------------------------------------------------------
# prepare data
#-----------------------------------------------------------------------------------
data(survival_swallows)
data <- survival_swallows
inits <- list(list(a=rep(3.5, len=data$K-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9))
data$ones <- matrix(1, nrow=data$I, ncol=data$K)
data$ones2 <- rep(1, len=data$I)
data$ind <- data$K:2                    # used to calculate chi
data$last <- as.vector(sapply(1:data$I, function(i) max(which(data$CH[i,]==1))))
saveRDS(data, 'data.RDS')
saveRDS(inits, 'inits.RDS')
params <- names(inits[[1]])
str(data)

## Run both models
.jags <- jags(data=data, model.file='swallows.jags',
              n.chains=1, n.iter=5000+5000, n.burnin=5000, n.thin=5,
              parameters.to.save=params)
.stan <- stan(file = "swallows.stan", data=data, chains=1, init=inits,
            iter=5000+5000, warmup=5000, pars=params, thin=5)
## get results and clean up for plotting
sims.jags <- as.data.frame(.jags$BUGSoutput$sims.matrix)
sims.stan <- as.data.frame(extract(.stan, permuted=FALSE)[,1,])
sims.stan$lp__ <- sims.jags$deviance <- NULL
## Look at effective sample sizes to make sure comparable
ess.jags <- monitor(sims=.jags$BUGSoutput$sims.array, warmup=0, print=FALSE)[,'n_eff'][params]
ess.stan <- monitor(sims=extract(.stan, permuted=FALSE), warmup=0, print=FALSE)[,'n_eff'][params]
data.frame(t(rbind(ess.stan, ess.jags)))

## Massage qqplot results into long format for ggplot
qq <- ldply(names(sims.jags), function(i){
                temp <- as.data.frame(qqplot(sims.jags[,i], sims.stan[,i], plot.it=FALSE))
                return(cbind(i,temp))})
g <- ggplot(qq, aes(x,y))+ geom_point(alpha=.5) +
    geom_abline(slope=1, col='red') + facet_wrap('i', scales='free') +
        xlab('jags')+ ylab('stan')
g
ggsave('model_comparison.png', g, width=9, height=7)
