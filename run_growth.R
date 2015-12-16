
Nfish.vec <- c(5, 10, 25, 50, 100, 200, 500)
L.vec <- c(1,5,20)
seeds <- c(4,6,1,90,15)
source("growth/run_chains.R")
source("growth_t/run_chains.R")

growth <- readRDS('growth/results/perf.RDS')
growth_t <- readRDS('growth_t/results/perf.RDS')
growth.df <- rbind(growth, growth_t)
growth.df$samples.per.time <- growth.df$minESS/growth.df$time
growth.df$minESS <- 100*growth.df$minESS/growth.df$Npar
growth.df.long <- melt(growth.df, c('model', 'platform', 'seed', 'Npar', 'Nsims'))
growth.df.long <- ddply(growth.df.long, .(platform, Npar, variable, model), mutate,
                       mean.value=mean(value))
growth.df.wide <- dcast(growth.df, platform+seed+Npar~model, value.var='time')
growth.df.wide$time.ratio <- growth.df.wide$growth_t/growth.df.wide$growth
## Use this to make function to compare between models later
g <- ggplot(growth.df.long, aes(Npar, log(value), color=platform)) +
    geom_jitter(position=position_jitter(width=.5, height=0), alpha=.5) +
        geom_line(data=growth.df.long, aes(Npar, log(mean.value))) +
            facet_grid(model~variable)
ggsave('plots/growth_comparison.png', g, width=ggwidth, height=ggheight)
g <- ggplot(growth.df.wide, aes(log(growth), log(growth_t), group=platform,
                           color=platform, size=platform=='jags')) +
    geom_point() + geom_abline(intercept=0, slope=1) + ggtitle("Time Comparison")
ggsave('plots/growth_comparison_time.png', g, width=ggwidth, height=ggheight)
g <- ggplot(growth.df.wide, aes(Npar, time.ratio, group=platform,
                           color=platform, size=platform=='jags')) +
    geom_point() + geom_abline(intercept=1, slope=0)+
        ggtitle("Time Comparison")
ggsave('plots/growth_comparison_timeratio.png', g, width=ggwidth, height=ggheight)
