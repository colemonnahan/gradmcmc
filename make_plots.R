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
g <- ggplot(adapt_empirical, aes(model, eps.final)) + geom_point(alpha=.5)
g <- ggplot(adapt_empirical, aes(model, delta.mean)) + geom_point(alpha=.5)
g <- ggplot(adapt_empirical, aes(model, nsteps.median)) + geom_point(alpha=.5)
g <- ggplot(adapt_empirical, aes(model, max_treedepths)) + geom_point(alpha=.5)
g <- ggplot(adapt_empirical, aes(eps.final, nsteps.median, color=model)) + geom_point()
