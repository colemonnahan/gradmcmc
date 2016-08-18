## Dig deeper into swallows model to see what is going on.

## Setup parallel
rstan_options(auto_write=TRUE)
ncores <- parallel::detectCores()-2
options(mc.cores=ncores)
### Manually read in the inputs for the two models
setwd('models/ss_logistic')
dat <- read.csv('tuna_data.csv')
cpue <- dat$CPUE
catches <- dat$Catches
data.ss_logistic<- list(catches=catches, logcpue=log(cpue), N=nrow(dat))
## sims.ind <- readRDS('sims.ind.RDS')
inits.ss_logistic.fn <- function()
  list(logr=log(.8), logK=log(279), iq=5, isigma2=100, itau2=100,
       u=rep(0, len=23))
pars.ss_logistic <- c('logr', 'logK', 'isigma2', 'iq', 'itau2', 'u')
setwd('../ss_logistic_nc/')
data.ss_logistic_nc <- data.ss_logistic
inits.ss_logistic_nc.fn <- function()
  list(logr=log(.8), logK=log(279), iq=5, isigma2=100, itau2=100,
       u_raw=rep(0, len=23))
## sims.ind <- readRDS('sims.ind.RDS')
pars.ss_logistic_nc <- c('logr', 'logK', 'isigma2', 'iq', 'itau2', 'u_raw')
setwd(main.dir)

### Run long chains with defaults to show divergences
stan.ss_logistic.fit <-
  stan(file='models/ss_logistic/ss_logistic.stan', par=pars.ss_logistic,
       init=inits.ss_logistic.fn, data=data.ss_logistic, iter=2000,
       chains=ncores, control=list(adapt_delta=.8))
my_sso <- shinystan::launch_shinystan(stan.ss_logistic.fit)
stan.ss_logistic_nc.fit <-
  stan(file='models/ss_logistic_nc/ss_logistic_nc.stan', par=pars.ss_logistic_nc,
       init=inits.ss_logistic_nc.fn, data=data.ss_logistic_nc, iter=2000,
       chains=ncores, control=list(adapt_delta=.8))
my_sso <- shinystan::launch_shinystan(stan.ss_logistic_nc.fit)
png('pairs_ss_logistic.png', res=500, width=7, height=5, units='in')
pairs(stan.ss_logistic.fit, pars=c('logK', 'logr'))
dev.off()
png('pairs_ss_logistic_nc.png', res=500, width=7, height=5, units='in')
pairs(stan.ss_logistic_nc.fit, pars=c('logK', 'logr'))
dev.off()

### Run highly thinned runs with STan and JAGS to try to get independent
### samples to test for bias
Nout.ind=200; Nthin.ind=2000
xx <- stan(file='ss_logistic.stan', par=params.jags, init=inits, data=data,
              iter=Nout.ind*Nthin.ind, thin=Nthin.ind, chains=10)
my_sso <- shinystan::launch_shinystan(xx)
temp <- extract(xx, permuted=FALSE)
dim(temp)<- c(dim(temp)[1]*dim(temp)[2], 1, dim(temp)[3])
perf.stan <- data.frame(rstan::monitor(sims=temp, warmup=0, print=FALSE, probs=.5))
sims.stan <- as.data.frame(temp);
names(sims.stan) <- names(xx)
yy <- jags.parallel(data=data, parameters.to.save=params.jags, inits=inits.fn,
           model.file='ss_logistic.jags', n.chains=10, n.iter=200*2000,
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
    stan.ss_logistic.fit <-
        stan(file='models/ss_logistic/ss_logistic.stan', par=pars.ss_logistic,
             init=inits.ss_logistic.fn, data=data.ss_logistic, iter=Niter,
             warmup=Nwarmup, chains=ncores, control=list(adapt_delta=ii))
    sims.stan <- extract(stan.ss_logistic.fit, permuted=FALSE)
    ## perf.stan.nuts <- monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5)
    adapt.nuts <- get_sampler_params(stan.ss_logistic.fit)
    for(i in 1:ncores){
        temp <- data.frame(adapt.nuts[[i]])
        adapt.list[[k]] <-
            data.frame(model='ss_logistic', chain=i,
                       delta.mean=mean(temp$accept_stat__[ind.samples]),
                       delta.target=ii, eps.final=tail(temp$stepsize__,1),
                       ndivergent=sum(temp$n_divergent__[ind.samples]))
        ## perf.list[[k]] <-
        ##   data.frame(model='ss_logistic', delta.target=ii,
        ##              eps.final=tail(adapt.nuts$stepsize__,1),
        ##              minESS=min(perf.stan.nuts$n_eff))
        k <- k+1
    }
    stan.ss_logistic_nc.fit <-
        stan(file='models/ss_logistic_nc/ss_logistic_nc.stan', par=pars.ss_logistic_nc,
             init=inits.ss_logistic_nc.fn, data=data.ss_logistic_nc, iter=Niter,
             warmup=Nwarmup, chains=ncores, control=list(adapt_delta=ii))
    sims.stan <- extract(stan.ss_logistic_nc.fit, permuted=FALSE)
    ## perf.stan.nuts <- monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5)
    adapt.nuts <- get_sampler_params(stan.ss_logistic_nc.fit)
    for(i in 1:ncores){
        temp <- data.frame(adapt.nuts[[i]])
        adapt.list[[k]] <-
            data.frame(model='ss_logistic_nc', chain=i,
                       delta.mean=mean(temp$accept_stat__[ind.samples]),
                       delta.target=ii, eps.final=tail(temp$stepsize__,1),
                       ndivergent=sum(temp$n_divergent__[ind.samples]))
        ## perf.list[[k]] <-
        ##   data.frame(model='ss_logistic_nc', delta.target=ii,
        ##              eps.final=tail(adapt.nuts$stepsize__,1),
        ##              minESS=min(perf.stan.nuts$n_eff))
        k <- k+1
    }
}
adapt <- do.call(rbind, adapt.list)
saveRDS(adapt, file='ss_logistic_tests/divergence_tests.RDS')
g <- ggplot(adapt, aes(delta.target, ndivergent, group=model, color=model)) +
    geom_point() + geom_line()
ggsave('ss_logistic_tests/divergence_tests.png', width=ggwidth, height=ggheight)
