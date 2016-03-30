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
verify <- TRUE
Nout <- 10000; Nthin <- 1; Nthin.ind <- 100
cor.vec <- c(0,1)
Npar.vec <- c(5, 50, 100, 200, 400, 1000)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
verify <- TRUE
Npar <- 5
Nout <- 500; Nthin <- 1; Nthin.ind <- 10
cor.vec <- c(0, .25, .5, .75, .8, .85, .9, .95, .99, .999, .9999)[c(1,11)]
Npar.vec <- c(2, 5)
source(paste0('models/',m,'/run_model.R'))

## Run growth tests, cross between centered/noncentered
Nout <- 500; Nthin <- 1; Nthin.ind <- 5
Npar.vec <- c(5,10,50, 100)[1:2]
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

## swallows; Example 14.5 from Korner-Nievergelt et al
Nout <- 10000; Nthin <- 1;
Nout.ind <- 200; Nthin.ind <- 100
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))

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
## empirical.means.normalized <-
##     ddply(empirical.means, .(platform, model), mutate,
##           normalized.samples.per.time=mean.samples.per.time/max(mean.samples.per.time))
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
                      mean.samples.per.time=mean(samples.per.time))
growth.means.wide <- dcast(subset(growth.means, Npar==15), model+Npar+centered+normal~platform, value.var='mean.samples.per.time')
growth.means.wide$stan_re <- with(growth.means.wide, round(stan.nuts/jags, 2))
mvn <- subset(simulated, model == 'mvn')
## Write them to file
write.csv(empirical, file='results/empirical.csv')
write.csv(simulated, file='results/simulated.csv')
### End of Step 3.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 4: Create plots, figures, and tables

m <- c('ss_logistic', 'redkite', 'growth_nct', 'swallows', 'mvnc','mnvd')
g <- ggplot(subset(empirical.means.normalized, platform!='jags' & model %in% m)) +
    geom_line(aes(delta.target, normalized.samples.per.time, color=model))
ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(all.wide, aes(log(Npar), log(stan.re.perf), color=model, shape=kind)) +
    geom_point() + geom_hline(yintercept=0)
ggsave('plots/perf_by_Npar.png', g, width=ggwidth, height=ggheight)
g <- ggplot(all.wide, aes(x=log(jags), log(stan.nuts), color=model, shape=kind,
                          size=Npar)) + geom_point() + geom_abline(1)
ggsave('plots/perf_by_Npar2.png', g, width=ggwidth, height=ggheight)
g <- ggplot(all, aes(log(minESS), y=log(minESS.coda), color=model))  +
    geom_point() + facet_wrap('platform') + xlim(0,7.5) +ylim(0,7.5) +
        geom_abline(intercept=0,slope=1)
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)
g <- ggplot(data=growth, aes(log(Npar), log(samples.per.time), color=platform)) +
    geom_point() + facet_grid(centered~normal)
ggsave('plots/perf_growth_simulated.png', g, width=ggwidth, height=ggheight)


g <- ggplot(mvn, aes(log(Npar), log(samples.per.time),
                          color=platform)) + geom_point()
ggsave('plots/perf_mvn_simulated.png', g, width=ggwidth,
       height=ggheight)

write.csv(file='results/table.csv', dcast(subset(empirical.means, platform=='jags' | delta.target==.8),
      formula=model~platform, value.var='mean.samples.per.time'))
write.csv(file='results/table_growth.csv', x=growth.means.wide)
### End of Step 4.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
