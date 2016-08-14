## Dig deeper into swallows model to see what is going on.

## Setup parallel
rstan_options(auto_write=TRUE)
ncores <- parallel::detectCores()-2
options(mc.cores=ncores)
### Manually read in the inputs for the two models
setwd('models/swallows')
data.swallows <- readRDS('data.RDS')
## sims.ind <- readRDS('sims.ind.RDS')
inits.swallows.fn <- function()
  list(a=rlnorm(17, log(3.5), .3), a1=runif(1, .2,1.5), b0=runif(4, -.5, 5),
       b1=runif(4,-.5, .1), sigmayearphi=runif(1, .01, 3),
       sigmaphi=runif(1,.01, 2), sigmap=runif(1, .3, 2), fameffphi=rep(0, len=72),
       fameffp=rep(0, len=72), yeareffphi=rep(0, len=4))
pars.swallows <-
    c('a', 'a1', 'b0', 'b1', 'sigmayearphi', 'sigmaphi', 'sigmap',
      'fameffphi', 'fameffp', 'yeareffphi')
setwd('../swallows_nc/')
data.swallows_nc <- readRDS('data.RDS')
## sims.ind <- readRDS('sims.ind.RDS')
inits.swallows_nc.fn <- function()
  list(a=rlnorm(17, log(3.5), .3), a1=runif(1, .2,1.5), b0=runif(4, -.5, 5),
       b1=runif(4,-.5, .1), sigmayearphi=runif(1, .01, 3),
       sigmaphi=runif(1,.01, 2), sigmap=runif(1, .3, 2), fameffphi_raw=rep(0, len=72),
       fameffp_raw=rep(0, len=72), yeareffphi_raw=rep(0, len=4))
pars.swallows_nc <-
    c('a', 'a1', 'b0', 'b1', 'sigmayearphi', 'sigmaphi', 'sigmap',
      'fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw')
setwd(main.dir)

### Run long chains with defaults to show divergences
stan.swallows.fit <-
  stan(file='models/swallows/swallows.stan', par=pars.swallows,
       init=inits.swallows.fn, data=data.swallows, iter=2000,
       chains=ncores, control=list(adapt_delta=.8))
my_sso <- shinystan::launch_shinystan(stan.swallows.fit)
stan.swallows_nc.fit <-
  stan(file='models/swallows_nc/swallows_nc.stan', par=pars.swallows_nc,
       init=inits.swallows_nc.fn, data=data.swallows_nc, iter=2000,
       chains=ncores, control=list(adapt_delta=.8))
my_sso <- shinystan::launch_shinystan(stan.swallows_nc.fit)
png('pairs_swallows.png', res=500, width=7, height=5, units='in')
pairs(stan.swallows.fit, pars=c('sigmayearphi', 'yeareffphi[1]'))
dev.off()
png('pairs_swallows_nc.png', res=500, width=7, height=5, units='in')
pairs(stan.swallows_nc.fit, pars=c('sigmayearphi', 'yeareffphi_raw[1]'))
dev.off()

### Run highly thinned runs with STan and JAGS to try to get independent
### samples to test for bias
Nout.ind=200; Nthin.ind=2000
xx <- stan(file='swallows.stan', par=params.jags, init=inits, data=data,
              iter=Nout.ind*Nthin.ind, thin=Nthin.ind, chains=10)
my_sso <- shinystan::launch_shinystan(xx)
temp <- extract(xx, permuted=FALSE)
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
perf.stan <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))
sims.stan <- as.data.frame(temp);
names(sims.stan) <- names(xx)
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


### Run Stan across a gradient of adapt_delta to see how %divergence changes
adapt.list <- perf.list <- list()
Niter <- 2*5000
Nwarmup <- Niter/2
ind.warmup <- 1:Nwarmup              # index of samples, excluding warmup
ind.samples <- (Nwarmup+1):Niter     # index of warmup samples
k <- 1
for(ii in c(.8,.85, .9, .95)) {
    stan.swallows.fit <-
        stan(file='models/swallows/swallows.stan', par=pars.swallows,
             init=inits.swallows.fn, data=data.swallows, iter=Niter,
             warmup=Nwarmup, chains=ncores, control=list(adapt_delta=ii))
    sims.stan <- extract(stan.swallows.fit, permuted=FALSE)
    ## perf.stan.nuts <- monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5)
    adapt.nuts <- get_sampler_params(stan.swallows.fit)
    for(i in 1:ncores){
        temp <- data.frame(adapt.nuts[[i]])
        adapt.list[[k]] <-
            data.frame(model='swallows', chain=i,
                       delta.mean=mean(temp$accept_stat__[ind.samples]),
                       delta.target=ii, eps.final=tail(temp$stepsize__,1),
                       ndivergent=sum(temp$n_divergent__[ind.samples]))
        ## perf.list[[k]] <-
        ##   data.frame(model='swallows', delta.target=ii,
        ##              eps.final=tail(adapt.nuts$stepsize__,1),
        ##              minESS=min(perf.stan.nuts$n_eff))
        k <- k+1
    }
    stan.swallows_nc.fit <-
        stan(file='models/swallows_nc/swallows_nc.stan', par=pars.swallows_nc,
             init=inits.swallows_nc.fn, data=data.swallows_nc, iter=Niter,
             warmup=Nwarmup, chains=ncores, control=list(adapt_delta=ii))
    sims.stan <- extract(stan.swallows_nc.fit, permuted=FALSE)
    ## perf.stan.nuts <- monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5)
    adapt.nuts <- get_sampler_params(stan.swallows_nc.fit)
    for(i in 1:ncores){
        temp <- data.frame(adapt.nuts[[i]])
        adapt.list[[k]] <-
            data.frame(model='swallows_nc', chain=i,
                       delta.mean=mean(temp$accept_stat__[ind.samples]),
                       delta.target=ii, eps.final=tail(temp$stepsize__,1),
                       ndivergent=sum(temp$n_divergent__[ind.samples]))
        ## perf.list[[k]] <-
        ##   data.frame(model='swallows_nc', delta.target=ii,
        ##              eps.final=tail(adapt.nuts$stepsize__,1),
        ##              minESS=min(perf.stan.nuts$n_eff))
        k <- k+1
    }
}
adapt <- do.call(rbind, adapt.list)
saveRDS(adapt, file='swallows_tests/divergence_tests.RDS')
g <- ggplot(adapt, aes(delta.target, ndivergent, group=model, color=model)) +
    geom_point() + geom_line()
ggsave('swallows_tests/divergence_tests.png', width=ggwidth, height=ggheight)
