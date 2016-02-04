### ------------------------------------------------------------
## This is the central file for running the MCMC efficiency analysis for my
## dissertation.
##
## Started 8/2015

### ------------------------------------------------------------
### Step 0: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'D:/gradmcmc/'
main.dir <- 'C:/Users/Cole/gradmcmc/'
setwd(main.dir)
source("startup.R")
Nout.ind <- 1000
seeds <- c(1:6)
lambda.vec <- NULL
delta.vec <- c(.5, .7, .8, .9, .95)
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
### End of Step 0.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 1

## Run multivariate normal, empirical and simulated
m <- 'mvn'
Nout <- 2000; Nthin <- 1; Nthin.ind <- 100
Npar.vec <- c(5, 10, 50, 100, 200, 500, 1000)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered and normal/t
## distributions
m <- 'growth'
Nout <- 10000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 100, 200, 500, 1000)
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_t'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nct'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

### End of Step 1.
### ------------------------------------------------------------
### ------------------------------------------------------------
### Step 2: Run simulations. For each model, run across the set of
### algorithms.
### End of Step 2.
### ------------------------------------------------------------



### ------------------------------------------------------------
### Step 3: Load and prepare data
perf.empirical <- ldply(list.files('results', pattern='perf_empirical'), function(i)
    read.csv(paste0('results/',i)))
perf.empirical.means <-
    ddply(perf.empirical, .(platform, model, delta.target), summarize,
          mean.samples.per.time=mean(samples.per.time))
perf.empirical.means <- ddply(perf.empirical.means, .(platform, model), mutate,
                              normalized.samples.per.time=mean.samples.per.time/max(mean.samples.per.time))
### End of Step 3.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 4: Create plots, figures, and tables
g <- ggplot(subset(perf.empirical.means, platform!='jags'),
            aes(delta.target, normalized.samples.per.time, color=model)) +
    geom_point() + geom_line()
ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)

write.csv(file='results/table.perf.csv', dcast(subset(perf.empirical.means, platform=='jags' | delta.target==.8),
      formula=model~platform, value.var='mean.samples.per.time'))
### End of Step 4.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
