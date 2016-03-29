## Load empirical data and inits
dat <- sample.lengths(Nfish=Npar, n.ages=5, logLinf.mean=logLinf.mean,
                           logLinf.sigma=logLinf.sigma, logk.mean=logk.mean,
                           logk.sigma=logk.sigma, sigma.obs=sigma.obs, t0=t0)
g <- ggplot(dat, aes(ages, loglengths, group=fish, color=fish)) +
    geom_point(alpha=.5) + geom_line()
ggsave('plots/simulated_growth.png', g, width=9, height=5)
data <- list(Nfish=Npar, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
inits <- list(logLinf_mean=logLinf.mean, logLinf_sigma=logLinf.sigma,
                  logk_mean=logk.mean, logk_sigma=logk.sigma, sigma_obs=sigma.obs,
                  logLinf=rep(logLinf.mean, len=Npar),
                  logk=rep(logk.mean, len=Npar), delta=1)
inits <- lapply(seq_along(seeds), function(i) inits)
