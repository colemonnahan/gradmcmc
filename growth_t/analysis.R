## Analyze the chains from run_chains.R

sims.jags.10 <- as.data.frame(readRDS('results/sims.jags.10.RDS')[1:50000,1,])
sims.jags.10.thinned <- sims.jags.10[seq(1, nrow(sims.jags.10), len=5000),]
sims.stan.nuts.10 <- as.data.frame(readRDS('results/sims.stan.nuts.10.RDS')[1:50000,1,])
sims.stan.nuts.10$lp__ <- NULL
sims.stan.nuts.10.thinned <- sims.stan.nuts.10[seq(1, nrow(sims.stan.nuts.10), len=5000),]
sims.stan.hmc1.10 <- as.data.frame(readRDS('results/sims.stan.hmc1.10.RDS')[1:50000,1,])
sims.stan.hmc1.10$lp__ <- NULL
sims.stan.hmc1.10.thinned <- sims.stan.hmc1.10[seq(1, nrow(sims.stan.hmc1.10), len=5000),]
sims.stan.hmc10.10 <- as.data.frame(readRDS('results/sims.stan.hmc10.10.RDS')[1:50000,1,])
sims.stan.hmc10.10$lp__ <- NULL
sims.stan.hmc10.10.thinned <- sims.stan.hmc10.10[seq(1, nrow(sims.stan.hmc10.10), len=5000),]


## make set of plots for each model to verify
make.trace <- function(df.thinned, model, string){
    nrows <- ceiling(sqrt(ncol(df.thinned)))
    png(paste0('plots/', model, '.trace.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df.thinned)){
        plot(df.thinned[,i], type='l', axes=FALSE,
             ylim=range(df.thinned[,i]), col=rgb(0,0,0,.5)); box()
        title(names(df.thinned)[i], line=-1)
    }
    dev.off()
}

make.trace(sims.jags.10.thinned, 'growth', 'jags.10')
make.acf(sims.jags.10, 'growth', 'jags.10')
make.trace(sims.stan.nuts.10.thinned, 'growth', 'stan.nuts.10')
make.acf(sims.stan.nuts.10, 'growth', 'stan.nuts.10')
make.trace(sims.stan.hmc1.10.thinned, 'growth', 'stan.hmc1.10')
make.acf(sims.stan.hmc1.10, 'growth', 'stan.hmc1.10')
make.trace(sims.stan.hmc10.10.thinned, 'growth', 'stan.hmc10.10')
make.acf(sims.stan.hmc10.10, 'growth', 'stan.hmc10.10')

make.acf <- function(df, model, string){
    nrows <- ceiling(sqrt(ncol(df)))
    png(paste0('plots/', model, '.acf.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df)) {
        acf(df[,i], axes=FALSE);box()
        title(names(df)[i], line=-1)
    }
    dev.off()
}



