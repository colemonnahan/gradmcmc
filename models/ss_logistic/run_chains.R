setwd('ss_logistic')
xx <- list()
for(Nyears in Nyears.vec){
    ## Reproducible data since seed set inside the function
    message(paste("======== Starting Nyears=", Nyears))
    set.seed(115)
    trash <- capture.output(source("load_models.R"))
    xx[[Nyears]] <-
        run.chains(model.name='ss_logistic', model.jags=model.jags, model.stan=model.stan,
                   inits.jags=inits.jags, params.jags=params.jags, data.jags=data.jags,
                   inits.stan=inits.stan, params.stan=params.stan, data.stan=data.stan,
                   seeds=seeds, Nout=Nout, n.burnin=n.burnin, n.thin=n.thin,
                   L=L.vec)
    ## Save them as we go in case it fails
    perf <- do.call(rbind, do.call(rbind, lapply(xx, function(x) x[['perf.list']])))
    trash <- do.call(rbind, lapply(xx, function(x) x[['adapt.list']]))
    adapt <- do.call(rbind.fill, trash[!sapply(trash, is.null)])
    saveRDS(perf, 'results/perf.RDS')
    saveRDS(adapt, 'results/adapt.RDS')
}
message("Making  plots")
perf <- readRDS('results/perf.RDS')
adapt <- readRDS('results/adapt.RDS')
plot.model.results(perf, adapt)
setwd('..')
