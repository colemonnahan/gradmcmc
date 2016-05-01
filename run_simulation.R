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
seeds <- c(1:10)
lambda.vec <- NULL
delta.vec <- .8 #c(.5, .7, .8, .9, .95)
metric <- c('unit_e', 'diag_e', 'dense_e')[2]
## Suppress JAGS and Stan output? Useful after debugging to clean up
## console and judge progress.
sink <- TRUE

.call('get_version', package='rjags')   # JAGS version 4.2.0
version$version.string                  # R version 3.2.3
packageVersion('rstan')                 # 2.8.2
packageVersion('R2jags')                # 0.5.7
packageVersion('rjags')                 # 4.4

### End of Step 0.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 1

## Run multivariate normal, empirical and simulated
m <- 'mvnd'
verify <- FALSE
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
## cor is a factor for independent (0) or from wishart (1)
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 25, 50, 100, 150, 200)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
verify <- FALSE
Npar <- 5
Nout <- 50000; Nthin <- 1; Nthin.ind <- 10
cor.vec <- c(0, .25, .5, .75, .8, .85, .9, .95, .99, .999)
Npar.vec <- c(2, 50, 100)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 100000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(5,10,50, 100)
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic
Nout <- 100000; Nthin <- 1; Nthin.ind <- 500
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 50000; Nthin <- 1; Nthin.ind <- 100
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows; Example 14.5 from Korner-Nievergelt et al
Nout <- 200000; Nthin <- 1; Nthin.ind <- 100
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

## quantgene; Example 14.5 from Korner-Nievergelt et al
Nout <- 2000; Nthin <- 1; Nthin.ind <- 100
m <- 'quantgene_nc'
source(paste0('models/',m,'/run_model.R'))
## m <- 'quantgene'
## source(paste0('models/',m,'/run_model.R'))


### End of Step 1.
### ------------------------------------------------------------
### ------------------------------------------------------------
### Step 2: Run simulations. For each model, run across the set of
### algorithms.
### End of Step 2.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables
setwd(main.dir)
## Loop through and get the maximum correlation of each model for merging
## into the main results.
cor.table <- ldply(list.files('models'), function(i) {
              xx <- readRDS(file.path('models', i, 'sims.ind.RDS'))
              cortemp <- cor(xx)
              max.cor <- max(abs(cortemp[lower.tri(cortemp)]))
              median.cor <- median(abs(cortemp[lower.tri(cortemp)]))
              data.frame(model=i, max.cor=max.cor, median.cor=median.cor)
            })
empirical <- ldply(list.files('results', pattern='perf_empirical'), function(i)
    read.csv(paste0('results/',i)))
## normalize by maximum run time across delta.target values
empirical <-
    ddply(empirical, .(platform, model, delta.target), mutate,
          mean.efficiency=mean(samples.per.time),
          sd.efficiency=sd(samples.per.time),
          median.efficiency=quantile(samples.per.time, probs=.5),
          lwr.efficiency=min(samples.per.time),
          upr.efficiency=max(samples.per.time))
empirical <- merge(x=empirical, y=cor.table, by='model')
empirical.means <-
    ddply(empirical, .(platform, model, Npar, max.cor), summarize,
          mean.samples.per.time=mean(samples.per.time))
empirical.means.wide <- dcast(subset(empirical.means), model+Npar+max.cor~platform,
                              value.var='mean.samples.per.time')
empirical.means.wide$stan_re <- with(empirical.means.wide, round(stan.nuts/jags, 2))
simulated <- ldply(list.files('results', pattern='perf_simulated'), function(i)
    read.csv(paste0('results/',i)))
simulated <-
    ddply(simulated, .(platform, model, delta.target, Npar, cor), mutate,
          mean.efficiency=mean(samples.per.time),
          sd.efficiency=sd(samples.per.time),
          median.efficiency=quantile(samples.per.time, probs=.5),
          lwr.efficiency=min(samples.per.time),
          upr.efficiency=max(samples.per.time))
## Select Stan modles with default delta.target level
growth <-
    subset(simulated, model %in%
               c('growth','growth_t','growth_nc','growth_nct'))
growth$model <- as.character(growth$model)
growth$centered <- 'centered'
growth$centered[growth$model %in% c('growth_nc', 'growth_nct')] <- 'noncentered'
growth$normal <- 'normal'
growth$normal[growth$model %in% c('growth_t', 'growth_nct')] <- 'student-t'
growth.means <- ddply(growth, .(platform, model, Npar, normal, centered), summarize,
                      mean.samples.per.time=round(mean(samples.per.time),2))
growth.means.wide <- dcast(subset(growth.means), model+Npar+centered+normal~platform, value.var='mean.samples.per.time')
growth.means.wide$stan_re <- with(growth.means.wide, round(stan.nuts/jags, 2))
mvn <- subset(simulated, model %in% c('mvnd', 'mvnc'))
mvn.means <- ddply(mvn, .(platform, model, Npar, cor), summarize,
                      mean.samples.per.time=round(mean(samples.per.time),2))
## Write them to file
write.table(empirical, file='results/empirical.csv', sep=',',
            row.names=FALSE, col.names=TRUE)
write.table(simulated, file='results/simulated.csv', sep=',',
            row.names=FALSE, col.names=TRUE)
write.table(file='results/table_growth.csv', x=growth.means.wide, sep=',',
            row.names=FALSE, col.names=TRUE)
write.table(file='results/table_cor.csv', x=cor.table, sep=',',
            row.names=FALSE, col.names=TRUE)
write.table(file='results/table_perf.csv', x=empirical.means.wide, sep=',',
            row.names=FALSE, col.names=TRUE)
print(subset(empirical, platform=='jags' & seed ==1, select=c(model, Nsims)))

### End of Step 3.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 4: Create exploratory plots
m <- c('ss_logistic', 'redkite', 'growth_nc', 'swallows', 'mvnc','mvnd')
## g <- ggplot(subset(empirical, platform!='jags' & model %in% m)) +
##     geom_point(aes(delta.target, log(samples.per.time))) + facet_wrap('model', scales='free_y')
## ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(minESS, y=(minESS-minESS.coda)/minESS))  +
    geom_point() + geom_abline(intercept=0,slope=0, col='red') +
      facet_grid(platform~model, scales='free')# + ylim(0,1)+ xlim(0,1)
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(model, y=(minESS/Nsims)))  +
  geom_point()+ ylim(0,.5) + facet_wrap('platform') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/ESS_percentages.png', g, width=ggwidth, height=ggheight)
g <- ggplot() + geom_line(data=growth.means,
                          aes(log10(Npar), log10(mean.samples.per.time),
                              group=centered, color=centered)) +
  geom_point(data=growth, aes(log10(Npar), y=log10(samples.per.time),
               color=centered)) + facet_grid(platform~normal)
ggsave('plots/perf_growth_simulated.png', g, width=ggwidth, height=ggheight)
g <- ggplot(subset(mvn.means, model=='mvnd'),
            aes(log10(Npar), log10(mean.samples.per.time), color=factor(cor), group=cor)) +
  geom_line() + facet_wrap('platform')
ggsave('plots/perf_mvnd_simulated.png', g, width=ggwidth,
       height=ggheight)
g <- ggplot(empirical.means.wide, aes(log10(Npar), max.cor, color=log10(stan_re)>0,
       size=abs(log10(stan_re)))) + geom_point()
ggsave('plots/perf_by_N_cor.png', g, width=ggwidth,
       height=ggheight)
### End of Step 4.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
test <- stan(file='swallows_nc.stan', data=data, init=inits,
             pars=params.jags, iter=2000, chains=1)

sims.stan.nuts <- extract(test, permuted=FALSE)
perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
names.temp <- row.names(perf.stan.nuts)[which(order(perf.stan.nuts$n_eff)<5)]

plot(test, pars=params.jags[1:7])
rstan::traceplot(test, pars=params.jags[1:7])
pairs(test, pars=names.temp)
library(shinystan)
launch_shinystan(test)

test <- stan(file='quantgene_nc.stan', data=data, init=inits,
             pars=pars, iter=2000, chains=1, control=list(adapt_delta=.9))
sims.stan.nuts <- extract(test, permuted=FALSE)
perf.stan.nuts <- data.frame(monitor(sims=sims.stan.nuts, warmup=0, print=FALSE, probs=.5))
eff <- sort(perf.stan.nuts$n_eff)
barplot(eff[1:50])
names.temp <- row.names(perf.stan.nuts)[which(order(perf.stan.nuts$n_eff)<15)]
plot(test, pars=pars[1:7])
rstan::traceplot(test, pars=names.temp)
pairs(test, pars=c('sa2',names.temp))
acf(sims.stan.nuts[,1, 'sm2'])
library(shinystan)
launch_shinystan(test)
