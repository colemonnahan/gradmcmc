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
m <- 'mvnd'
Nout <- 2000; Nthin <- 1; Nthin.ind <- 100
Npar.vec <- c(5, 10, 50, 100, 200, 500, 1000)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
Npar <- 5
Nout <- 10000; Nthin <- 1; Nthin.ind <- 10
cor.vec <- seq(0,.99, len=5)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 10000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 100, 200, 500, 1000)
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_t'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nct'
source(paste0('models/',m,'/run_model.R'))

## State space logistic
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 20000; Nthin <- 1; Nthin.ind <- 1000
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## Swiss tit example from Kery and Schaub; 4 full model w/ informative
## prior. BROKEN BECAUSE THE DATA INPUTS ARE DIFFERENT DUE TO MISSING VALUES
Nout <- 2000; Nthin <- 1; Nthin.ind <- 10
m <- 'tits'
source(paste0('models/',m,'/run_model.R'))

## Example 14.5 from Korner-Nievergelt et al
Nout <- 5000; Nthin <- 1; Nthin.ind <- 1
m <- 'swallows'
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
## normalize by maximum run time across delta.target values
perf.empirical.means <-
    ddply(perf.empirical, .(platform, model, delta.target), mutate,
          median.samples.per.time=median(samples.per.time))
perf.empirical.means <-
    ddply(perf.empirical.means, .(platform, model), mutate,
          normalized.samples.per.time=samples.per.time/max(median.samples.per.time))
perf.empirical.means <-
    ddply(perf.empirical.means, .(platform, model, delta.target), mutate,
          med=quantile(normalized.samples.per.time, probs=c(.5)),
          upr=quantile(normalized.samples.per.time, probs=c(.75)),
          lwr=quantile(normalized.samples.per.time, probs=c(.25)))

perf.simulated <- ldply(list.files('results', pattern='perf_simulated'), function(i)
    read.csv(paste0('results/',i)))
perf.all <- rbind(cbind(perf.empirical, kind='empirical'),
                  cbind(perf.simulated, kind='simulated'))
perf.all.wide <-
    dcast(subset(perf.all, platform=='jags' | delta.target==.8),
          kind+Npar+seed+model~platform, value.var='samples.per.time')
perf.all.wide <- within(perf.all.wide, stan.re.perf<-stan.nuts/jags)
perf.growth <-
    subset(perf.simulated, model %in% c('growth','growth_t','growth_nc','growth_nct'))
perf.growth$model <- as.character(perf.growth$model)
perf.growth$centered <- 'centered'
perf.growth$centered[perf.growth$model %in% c('growth_nc', 'growth_nct')] <- 'noncentered'
perf.growth$normal <- 'normal'
perf.growth$normal[perf.growth$model %in% c('growth_t', 'growth_nct')] <- 'student-t'
perf.mvn <- subset(perf.simulated, model == 'mvn')

### End of Step 3.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 4: Create plots, figures, and tables

m <- c('ss_logistic', 'redkite')
g <- ggplot(subset(perf.empirical.means, platform!='jags' & model %in% m)) +
 geom_point(aes(delta.target, normalized.samples.per.time, color=model))+
    geom_linerange(aes(delta.target, ymax=upr, ymin=lwr, color=model)) +
    geom_line(aes(delta.target, med, color=model))
g


ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(perf.all.wide, aes(Npar, log(stan.re.perf), color=model, shape=kind)) +
    geom_point() + geom_hline(yintercept=0)
ggsave('plots/perf_by_Npar.png', g, width=ggwidth, height=ggheight)
g <- ggplot(perf.all.wide, aes(log(jags), log(stan.nuts), color=model, shape=kind,
                          size=Npar)) + geom_point() + geom_abline(1)
ggsave('plots/perf_by_Npar2.png', g, width=ggwidth, height=ggheight)
g <- ggplot(perf.growth, aes(log(Npar), log(samples.per.time),
                             color=platform)) + geom_point() +
                                 facet_grid(centered~normal)
ggsave('plots/perf_growth_simulated.png', g, width=ggwidth, height=ggheight)
g <- ggplot(perf.mvn, aes(log(Npar), log(samples.per.time),
                          color=platform)) + geom_point()
ggsave('plots/perf_mvn_simulated.png', g, width=ggwidth,
       height=ggheight)

write.csv(file='results/table.perf.csv', dcast(subset(perf.empirical.means, platform=='jags' | delta.target==.8),
      formula=model~platform, value.var='mean.samples.per.time'))
### End of Step 4.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
