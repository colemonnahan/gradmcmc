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
Nout <- 100000; Nthin <- 1; Nthin.ind <- 100
Npar.vec <- c(5,10,50, 100)
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_t'
source(paste0('models/',m,'/run_model.R'))
m <- 'growth_nct'
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
print(subset(empirical, platform=='jags' & seed ==1, select=c(model, Nsims)))
### End of Step 3.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 4: Create exploratory plots
|m <- c('ss_logistic', 'redkite', 'growth_nct', 'swallows', 'mvnc','mnvd')
g <- ggplot(subset(empirical, platform!='jags' & model %in% m)) +
    geom_point(aes(delta.target, log(samples.per.time))) + facet_wrap('model', scales='free_y')
ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes((minESS/Nsims), y=(minESS.coda/Nsims), color=model))  +
    geom_point() + geom_abline(intercept=0,slope=1) +
      facet_grid(platform~model) + ylim(0,1)+ xlim(0,1)
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(model, y=(minESS/Nsims)))  +
  geom_violin()+ ylim(0,1) + facet_wrap('platform') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/ESS_percentages.png', g, width=ggwidth, height=ggheight)
g <- ggplot(data=growth.means,
            aes(log(Npar), log(mean.samples.per.time), group=centered,
                color=centered)) +
    geom_line(alpha=.5) + facet_grid(platform~normal)
ggsave('plots/perf_growth_simulated.png', g, width=ggwidth, height=ggheight)
g <- ggplot(subset(mvn.means, model=='mvnd'),
            aes(log(Npar), log(mean.samples.per.time), color=factor(cor), group=cor)) +
  geom_line() + facet_wrap('platform')
ggsave('plots/perf_mvn_simulated.png', g, width=ggwidth,
       height=ggheight)
### End of Step 4.
### ------------------------------------------------------------



### ------------------------------------------------------------
## Development code
