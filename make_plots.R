m <- c('ss_logistic', 'wildflower_nc', 'redkite', 'growth_nc', 'swallows', 'mvnc','mvnd')
## g <- ggplot(subset(empirical, platform!='jags' & model %in% m)) +
##     geom_point(aes(delta.target, log(samples.per.time))) + facet_wrap('model', scales='free_y')
## ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)


g <- ggplot(empirical, aes(minESS, y=minESS.coda, group=platform, color=platform))  +
  scale_x_log10()+scale_y_log10()+
    geom_point() + geom_abline(intercept=0,slope=1)+
      facet_wrap('model', scales='free')
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)

g <- ggplot(empirical, aes(model, y=(minESS/Nsims)))  +
  geom_point()+ ylim(0,.5) + facet_wrap('platform') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/ESS_percentages.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(model, y=log(time)))  +
  geom_point()+ facet_wrap('platform') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/runtime.png', g, width=ggwidth, height=ggheight)
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
empirical.means.wide$stan_faster <- empirical.means.wide$stan_re>1
empirical.means.wide$stan_faster <- empirical.means.wide$stan_re>1
g <- ggplot(empirical.means.wide, aes(log10(Npar), max.cor)) +
  geom_point(aes(size=abs(log10(stan_re)), color=stan_faster)) +
    geom_text(aes(label=model), hjust=1, vjust=2)
ggsave('plots/perf_by_N_cor.png', g, width=ggwidth,
       height=ggheight)
## Some adaptation plots
g <- ggplot(adapt_empirical, aes(model, eps.final)) + geom_point(alpha=.5)
ggsave('plots/adapt_eps.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, delta.mean)) + geom_point(alpha=.5)
ggsave('plots/adapt_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, log10(nsteps.median))) + geom_point(alpha=.5)
ggsave('plots/adapt_nsteps.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, max_treedepths)) + geom_point(alpha=.5)
ggsave('plots/adapt_max_treedepths.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, 100*ndivergent/Nsims)) +
    geom_point(alpha=.5) + ylab('Percent diverged transitions')
ggsave('plots/adapt_ndivergent.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(eps.final, log10(nsteps.mean), color=model)) + geom_point()
ggsave('plots/adapt_eps_vs_nsteps.png', g, width=ggwidth, height=ggheight)

ggplot(empirical, aes(time.warmup, time.sampling, color=model)) +
  facet_wrap('platform') + geom_point()
