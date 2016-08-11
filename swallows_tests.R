## Dig deeper into swallows model to see what is going on.

## Setup parallel
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
Nout.ind=200; Nthin.ind=2000
setwd('models/swallows')
data <- readRDS('data.RDS')
inits <- list(list(a=rep(3.5, len=data$K-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9,
                   fameffphi=rep(0, len=data$nfam),
                   fameffp=rep(0, len=data$nfam),
                   yeareffphi=rep(0, len=4)))
inits <- rep(inits, times=10)
## JAGS complains about environments so easiest is to pass everything as
## absolute numbers.
inits.fn <- function() list(a=rep(3.5, len=18-1), a1=0, b0=rep(2, len=4), b1=rep(0, len=4),
                   sigmayearphi=.7, sigmaphi=.5, sigmap=.9,
                   fameffphi=rep(0, len=72),
                   fameffp=rep(0, len=72),
                   yeareffphi=rep(0, len=4))
params.jags <-
    c('a', 'a1', 'b0', 'b1', 'sigmayearphi', 'sigmaphi', 'sigmap',
      'fameffphi', 'fameffp', 'yeareffphi')

stan.swallows.fit <- stan(file='swallows.stan', par=params.jags, init=inits, data=data,
              iter=20000, thin=10, chains=10, control=list(adapt_delta=.95))
my_sso <- shinystan::launch_shinystan(stan.fit)

xx <- stan(file='swallows.stan', par=params.jags, init=inits, data=data,
              iter=Nout.ind*Nthin.ind, thin=Nthin.ind, chains=10)
my_sso <- shinystan::launch_shinystan(xx)
temp <- extract(xx, permuted=FALSE)
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
perf.stan <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))
sims.stan <- as.data.frame(temp);
names(sims.stan) <- names(xx)

list2env(data, envir=globalenv())
list2env(inits, envir=globalenv())
yy <- jags.parallel(data=data, parameters.to.save=params.jags, inits=inits.fn,
           model.file='swallows.jags', n.chains=10, n.iter=200*2000,
           n.thin=2000)
temp <- yy$BUGSoutput$sims.array
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
sims.jags <- as.data.frame(yy$BUGSoutput$sims.matrix)
perf.jags <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))

perf.platforms <- rbind(cbind(platform='jags',perf.jags),
                        cbind(platform='stan',perf.stan))
perf.platforms <- melt(perf.platforms, c('Rhat', 'n_eff'), id.vars='platform')
plot.model.comparisons(sims.stan, sims.jags, perf.platforms)
