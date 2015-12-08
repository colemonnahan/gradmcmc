## ehook: An analysis of the Hamley and Skud (1978) data on the effect of
## hook

setwd('growth')
xx <- list()
file.remove('sink_progress.txt')
for(Nfish in Nfish.vec){
    ## Reproducible data since seed set inside the function
    message(paste("======== Starting Nfish=", Nfish))
    set.seed(115)
    source("load_models.R")
    xx[[Nfish]] <-
        run.chains(model.name='growth', model.jags=model.jags, model.stan=model.stan,
                   inits.jags=inits.jags, params.jags=params.jags, data.jags=data.jags,
                   inits.stan=inits.stan, params.stan=params.stan, data.stan=data.stan,
                   seeds=seeds, Nout=Nout, n.burnin=n.burnin, n.thin=n.thin,
                   L=L.vec)
    ## Save them as we go in case it fails
    perf.list <- lapply(xx, function(x) x[['perf.list']])
    adapt.list <- lapply(xx, function(x) x[['adapt.list']])
    saveRDS(perf.list, 'results/perf.list.RDS')
    saveRDS(adapt.list, 'results/adapt.list.RDS')
}

message("Making  plots")
perf.list <- readRDS('results/perf.list.RDS')
adapt.list <- readRDS('results/adapt.list.RDS')
plot.model.results(perf.list, adapt.list)

setwd('..')





