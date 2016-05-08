## Sourcing this file will run everything for this model, given the MCMC
## arguments are in the global workspace.
setwd(paste0('models/',m))

## Load empirical data and inits
dat <- read.csv('tuna_data.csv')
cpue <- dat$CPUE
catches <- dat$Catches
data <- list(catches=catches, logcpue=log(cpue), N=nrow(dat))
inits <- list(list(logr=log(.8), logK=log(279), iq=5, isigma2=100, itau2=100,
             u=rep(0, len=nrow(dat))))
params.jags <- c('logr', 'logK', 'isigma2', 'iq', 'itau2', 'u')
r <- exp((inits[[1]]$logr))
K <- exp((inits[[1]]$logK))
sd.process <- sqrt(1/inits[[1]]$isigma2)
trajectory <- rep(NA, data$N)
trajectory[1] <- K*exp(inits[[1]]$u[1]-sd.process^2/2)
for( yr in 2: data$N){
    Nt <- trajectory[yr-1]
    Ct <- data$catches[yr-1]
    trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct)*exp(inits[[1]]$u[yr]-sd.process^2/2)
}
png('plots/initial_fits.png', width=7, height=5, res=300, units='in')
plot(trajectory, type='l', ylim=c(0, 1.5*max(trajectory)))
points(exp(data$logcpue)*inits[[1]]$iq, col='red')
dev.off()

## Get independent samples from each model to make sure they are coded the
## same
if(verify)
verify.models(model=m, params.jags=params.jags, inits=inits, data=data,
              Nout=Nout.ind, Nthin=Nthin.ind)

sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
inits <- lapply(1:length(seeds), function(i)
  list(logr=sims.ind$logr[i],
       logK=sims.ind$logK[i],
       isigma2=sims.ind$isigma2[i],
       itau2=sims.ind$itau2[i],
       iq=sims.ind$iq[i],
        u=as.numeric(sims.ind[i, grep('u\\.', x=names(sims.ind))])))

## Fit empirical data with no thinning for efficiency tests
fit.empirical(model=m, params.jag=params.jags, inits=inits, data=data,
              lambda=lambda.vec, delta=delta, metric=metric, seeds=seeds,
              Nout=Nout)

### NOT USED!
## ## Now loop through model sizes and run for default parameters, using JAGS
## ## and NUTS only.
## adapt.list <- perf.list <- list()
## for(i in seq_along(Npar.vec)){
##     Npar <- Npar.vec[i]
##     ## Reproducible data since seed set inside the function
##     message(paste("======== Starting Npar=", Npar))
##     set.seed(115)
##     source("generate_data.R")
##     temp <- run.chains(model=m, inits=inits, params.jags=params.jags, data=data,
##                    seeds=seeds, Nout=Nout, Nthin=1, lambda=NULL)
##     adapt.list[[i]] <- temp$adapt
##     perf.list[[i]] <- temp$perf
##     ## Save them as we go in case it fails
##     perf <- do.call(rbind, perf.list)
##     adapt <- do.call(rbind, adapt.list)
##     write.csv(x=perf, file=results.file(paste0(m,'_perf_simulated.csv')))
##     write.csv(x=adapt, file=results.file(paste0(m,'_adapt_simulated.csv')))
##     rm(temp)
## }
## message("Making plots...")
## plot.simulated.results(perf, adapt)

setwd(main.dir)
