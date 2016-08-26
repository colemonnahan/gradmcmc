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
print(subset(empirical, platform=='jags' & seed ==seeds[1],
             select=c(model, Nsims)))
print(ddply(empirical, .(model, platform), summarize, reps=length(seeds)))

## summarize adaptation info
adapt_empirical<- ldply(list.files('results', pattern='adapt_empirical'), function(i)
    read.csv(paste0('results/',i)))
