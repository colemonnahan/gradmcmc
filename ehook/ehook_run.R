## An analysis of the Hamley and Skud (1978) data on the effect of hook
## spacing.

## The original data, to be processed into the correct format for each
## model below
setwd('ehook')
hs <- read.csv('data/hs_data.csv')
spacing <- hs$spacing
group <- as.numeric(hs$group)
## The dependent variable (catch.per.hook for now) on log scale since it
## needs to be positive
yobs <- hs$catch.per.hook
log.yobs <- log(yobs+.1)
Nobs <- length(yobs)
Ngroup <- length(unique(group))
day <- hs$daynumber

## Got some reasonable values from an initial run
inits <- list(beta=.1, gamma=.05, logcpue=rep(1,len=Ngroup),
              logcpue_mean=.2, logcpue_sd=.5, logsigma_obs=rep(-.5, len=Ngroup),
              sigma_obs_mean=-.5, sigma_obs_sd=.3)


## These are global parameters. Each software platform specifies them
## slightly different so see the function calls.
n.out1 <- 100
n.thin1 <- 2
n.chains1 <- 1
n.iter1 <- 1.25*n.out1*n.thin1
n.burnin1 <- .2*n.iter1

n.out.ind <- 1000
n.thin.ind <- 1000
n.chains.ind <- 4
n.iter.ind <- 1.25*n.out.ind*n.thin.ind
n.burnin.ind <- .2*n.iter.ind

### ------------------------------------------------------------
## Run the JAGS models
data.jags <-
    list(log_yobs=log.yobs, group=group, Nobs=Nobs, Ngroup=Ngroup, day=day,
         spacing=spacing)
params.jags <-
    c("logcpue_mean", "logcpue_sd", "logsigma_obs", "gamma", "logcpue",
      "sigma_obs_mean", "sigma_obs_sd", "beta")
## JAGS models
model.jags <- function(){
### The Hamley and Skud formula with a random effect on CPUE. Here there is
### a random effect on the observation error
    ## hyperpriors
    logcpue_mean~dunif(-5,5)            # Mean CPUE (asymptote)
    logcpue_sd~dunif(0,5)             # SD of mean CPUE
    sigma_obs_mean~dunif(-1,5)               # mean Observation error, log scale
    sigma_obs_sd~dunif(0,3)               # sd of obs errors, log scale
    beta~dunif(0,10)
    ## priors
    gamma~dunif(0,1)                   # Impact of day on CPUE
    ## Loop through the hyperparameters (on group) and calculate
    ## probabilities
    for(i in 1:Ngroup){
        logcpue[i]~dnorm(logcpue_mean, pow(logcpue_sd, -2))
        logsigma_obs[i]~dnorm(sigma_obs_mean, pow(sigma_obs_sd, -2))
    }
    ## Loop through observations and calculate likelihood
    for(i in 1:Nobs){
        ypred[i] <- exp(logcpue[group[i]])*exp(-day[i]*gamma)*(1-exp(-beta*spacing[i]))
        ## Likelihood of data
        log_yobs[i]~dnorm(log(ypred[i]), pow(exp(logsigma_obs[group[i]]), -2))
    }
}

## Run a long chain with thinning to get independent samples to make sure
## the models are matching.
results.jags.ind <- jags(data=data.jags, parameters.to.save=params.jags,
                           model.file=model.jags, n.chains=n.chains.ind, n.iter=n.iter.ind,
                           n.burnin=n.burnin.ind, n.thin=n.thin.ind)
saveRDS(results.jags.ind, file='results/results.jags.ind.RDS')
results.jags <- jags(data=data.jags, parameters.to.save=params.to.save.jags,
                           model.file=mod.jags, n.chains=n.chains1, n.iter=n.iter1,
                           n.burnin=n.burnin1, n.thin=n.thin1)
saveRDS(results1.jags, file='results/results1.jags.RDS')
## End of JAGS runs for this model
### ------------------------------------------------------------

### ------------------------------------------------------------
## Stan runs
data.stan <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
## Run a dummy chain to get the compilation out of the way for more
## sensible speed comparisons
##if(!exists('stan.model'))
stan.model <- stan(file='ehook.stan', data=data.stan, iter=1, chains=1, init=list(inits))
results.stan.ind <-
    stan(fit=stan.model, data=data.stan, iter=500000, warmup=2000,
         chains=5, thin=10, algorithm='NUTS', init=rep(list(inits),5),
         control=list(adapt_engaged=TRUE))
traceplot(results.stan.ind, pars=c("logcpue_mean", "logcpue_sd", "beta",
                            "gamma", "lp__"))
plot(results.stan.ind)
xx <- data.frame(extract(results.stan.ind))
saveRDS(results.stan.ind, file='results/results.stan.ind.RDS')


### ------------------------------------------------------------
## devtools::install_github('kaskr/adcomp/TMB')
library(TMB)
## TMB runs
data.tmb <-
    list(Ngroup=Ngroup, Nobs=Nobs,log_yobs=log.yobs, group=group, day=day,
         spacing=spacing)
data.tmb$group <- data.tmb$group-1
## Need to massage inits b/c of transformed parameters
inits.tmb <- inits
boundpinv <- function(x, min, max){
    -log( (max-min)/(x-min) -1)
}
inits.tmb$logcpue_mean2=boundpinv(inits.tmb$logcpue_mean, -5, 5)
inits.tmb$logcpue_sd2=boundpinv(inits.tmb$logcpue_sd, 0, 5)
inits.tmb$sigma_obs_mean2=boundpinv(inits.tmb$sigma_obs_mean, -5, 5)
inits.tmb$sigma_obs_sd2=boundpinv(inits.tmb$sigma_obs_sd, 0,5)
inits.tmb$beta2= boundpinv(inits.tmb$beta, 0,10)
inits.tmb$gamma2=boundpinv(inits.tmb$gamma, 0,1)
inits.tmb$logcpue_mean=NULL
inits.tmb$logcpue_sd=NULL
inits.tmb$sigma_obs_mean=NULL
inits.tmb$sigma_obs_sd=NULL
inits.tmb$beta=NULL
inits.tmb$gamma=NULL

## Need to reorder parameter list for TMB.
pars <- c("beta2","gamma2","logcpue", "logcpue_mean2","logcpue_sd2",
          "logsigma_obs", "sigma_obs_mean2", "sigma_obs_sd2" )
inits.tmb <- inits.tmb[pars]
compile("ehook.cpp")
dyn.load(TMB::dynlib("ehook"))
model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='ehook')
model.tmb$hessian <- TRUE
model.tmb.opt <- do.call(optim, model.tmb)
est.tmb <- model.tmb.opt$par
covar.tmb <- solve(model.tmb.opt$hessian)




### quick check that the adaption is going well and long enough
## test <- llply(delta.vec, function(x)
##     run_nuts(obj=model.tmb, nsim=nsim, seed=seed, covar=cov,
##              Madapt=Madapt, delta=x, diag=TRUE, inits=inits.tmb))
## test2 <- data.frame(do.call(cbind, llply(test, function(x) x$epsbar)))
## names(test2) <- delta.vec
## test2$iteration <- 1:nrow(test2)
## ggplot(reshape2::melt(test2, c('iteration')),
##        aes(iteration, value, group=variable, color=variable)) +
##     geom_line()

library(snowfall)
Ncores <- 5
sfStop()
sfInit( parallel=TRUE, cpus=Ncores )
thin <- 1
warmup <- 200
Nind <- 2000
Nper <- ceiling(Nind/Ncores)
Nout <- (Nper*thin+warmup)
sfExport("data.tmb", "inits.tmb", "Nout", "warmup", "est.tmb")
temp <- sfLapply(1:Ncores, function(i){
  dyn.load(TMB::dynlib('ehook'))
  model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='ehook')
  set.seed(i)
  x <- TMB::mcmc(obj=model.tmb, nsim=Nout, max_doubling=7,
            Madapt=warmup, delta=.5, algorithm='NUTS', params.init=est.tmb)
  nll <- apply(x, 1, model.tmb$fn)
  cbind(seed=i, x, nll)
})
## discard warmup and then thin
temp2 <- lapply(temp, function(x) x[-(1:warmup),])
temp3 <- lapply(temp2, function(x) x[seq(1, nrow(x), by=thin),])
nlls <- do.call(rbind, lapply(1:Ncores, function(x) data.frame(iteration=1:nrow(temp3[[x]]),seed=x, nll=temp3[[x]]$nll )))
ggplot(nlls, aes(iteration, nll, group=seed, color=factor(seed)))+geom_line()
tmb.ind.ess <- data.frame(do.call(rbind, lapply(temp3, effectiveSize)))
tmb.ind.ess$seed <- 1:Ncores
tmb.ind.ess.long <- reshape2::melt(tmb.ind.ess, 'seed', value.name='ess')
ggplot(tmb.ind.ess.long, aes(variable, ess, color=factor(seed)))+geom_point()

results.tmb.ind <- do.call(rbind, temp3)
## results.tmb.ind$LP <- -apply(results.tmb.ind[,-1], 1, model.tmb$fn)
saveRDS(results.tmb.ind, file='results/results.tmb.ind.RDS')

delta.vec <- seq(.1,.9, len=Ncores*2)
seeds.vec <- 1
nsim <- 10000
Madapt <- min(200, .1*nsim)
params <- data.frame(expand.grid(delta=delta.vec, covar=c(TRUE, FALSE), seed=seeds.vec))
sfExport("params", "Madapt", "nsim", "covar.tmb", "run_nuts", "data.tmb",
         "inits.tmb", "run_hmc")
ehook.nuts.perf.list <- sfLapply(1:nrow(params), function(i){
  dyn.load(TMB::dynlib("ehook"))
  model.tmb <- TMB::MakeADFun(data.tmb, parameters=inits.tmb, DLL='ehook')
   covar.temp <- if(params$covar[i]) covar.tmb else NULL
   run_nuts(obj=model.tmb, nsim=nsim, seed=params$seed[i],
            covar=covar.temp, Madapt=Madapt, delta=params$delta[i])
})
ehook.nuts.perf <- do.call(rbind, ehook.nuts.perf.list)
saveRDS(ehook.nuts.perf, results.file('ehook.nuts.perf.RDS'))
ehook.nuts.perf$acceptance <- NULL
ehook.nuts.perf.long <-
    reshape2::melt(ehook.nuts.perf, c('covar', 'tuning', 'algorithm', 'seed'))
ggplot(ehook.nuts.perf.long, aes(tuning, value, group=covar, color=covar)) +
    geom_line() + facet_wrap('variable', nrow=3, scales='free')
ggsave(plots.file('ehook_nuts_perf.png'), width=8, height=6)



### Old way of using plyr, but can't get parallel to work
## ehook.nuts.perf <-
##     ldply(seeds.vec, function(seed)
##           ldply(covar.list, function(cov)
##               llply(delta.vec, function(x)
##                   run_nuts(obj=model.tmb, nsim=nsim, seed=seed, covar=cov,
##                            Madapt=Madapt, delta=x)
##                     )))





function (posterior, mle, diag = c("acf", "hist", "trace"), acf.ylim = c(-1,
    1), ymult = NULL, axis.col = gray(0.5), which.keep = NULL,
    label.cex = 0.5, limits = NULL, ...)
{
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    diag <- match.arg(diag)
    if (NCOL(posterior) != mle$nopar)
        stop("Number of parameters in posterior and mle not the same")
    if (is.null(which.keep))
        which.keep <- 1:NCOL(posterior)
    par.names <- mle$names[which.keep]
    n <- length(par.names)
    if (n == 1)
        stop("This function is only meaningful for >1 parameter")
    mle.par <- mle$est[which.keep]
    mle.se <- mle$std[which.keep]
    mle.cor <- mle$cor[which.keep, which.keep]
    posterior <- posterior[, which.keep]
    if (is.null(ymult))
        ymult <- rep(1.3, n)
    if (is.null(limits)) {
        limits <- list()
        for (i in 1:n) {
            limit.temp <- mle.par[i] + c(-1, 1) * 1.96 * mle.se[i]
            min.temp <- min(posterior[, i], limit.temp[1])
            max.temp <- max(posterior[, i], limit.temp[2])
            margin <- 0.15 * (max.temp - min.temp)
            limits[[i]] <- c(min.temp - margin, max.temp + margin)
        }
    }
    par(mfrow = c(n, n), mar = 0 * c(0.1, 0.1, 0.1, 0.1), yaxs = "i",
        xaxs = "i", mgp = c(0.25, 0.25, 0), tck = -0.02, cex.axis = 0.65,
        col.axis = axis.col, oma = c(2, 2, 2, 0.5))
    temp.box <- function() box(col = axis.col, lwd = 0.5)
    for (row in 1:n) {
        for (col in 1:n) {
            if (row == col) {
                if (diag == "hist") {
                  h <- hist(posterior[, row], plot = F)
                  if (is.null(limits)) {
                    hist(posterior[, row], axes = F, freq = FALSE,
                      ann = F, ylim = c(0, ymult[row] * max(h$density)),
                      col = gray(0.8), border = gray(0.5))
                  }
                  else {
                    hist(posterior[, row], axes = F, freq = FALSE,
                      ann = F, ylim = c(0, ymult[row] * max(h$density)),
                      col = gray(0.8), border = gray(0.5), xlim = limits[[row]])
                  }
                  temp.box()
                }
                else if (diag == "acf") {
                  acf(posterior[, row], axes = F, ann = F, ylim = acf.ylim)
                  temp.box()
                }
                else if (diag == "trace") {
                  plot(x = posterior[, row], lwd = 0.5, col = gray(0.5),
                    type = "l", axes = F, ann = F, ylim = limits[[row]])
                  temp.box()
                }
                mtext(par.names[row], line = -2, cex = label.cex)
            }
            if (row > col) {
                par(xaxs = "r", yaxs = "r")
                plot(x = posterior[, col], y = posterior[, row],
                  axes = FALSE, ann = FALSE, pch = ifelse(NROW(posterior) >=
                    5000, ".", 1), xlim = limits[[col]], ylim = limits[[row]],
                  ...)
                points(x = mle.par[col], y = mle.par[row], pch = 16,
                  cex = 0.1, col = 2)
                ellipse.temp <- ellipse::ellipse(x = mle.cor[col,
                  row], scale = mle.se[c(col, row)], centre = mle.par[c(col,
                  row)], npoints = 1000, level = 0.95)
                lines(ellipse.temp, lwd = 1.5, lty = 1, col = "red")
                par(xaxs = "i", yaxs = "i")
                temp.box()
            }
            if (row < col) {
                plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0,
                  1), axes = F, ann = F)
                temp.cor <- round(cor(posterior[, c(row, col)])[1,
                  2], 2)
                legend("center", legend = NA, title = temp.cor,
                  cex = (3 * abs(temp.cor) + 0.25) * 0.5, bty = "n")
                temp.box()
            }
            if (row == n) {
                par(mgp = c(0.05, ifelse(col%%2 == 0, 0, 0.5),
                  0))
                axis(1, col = axis.col, lwd = 0.5)
            }
            if (col == 1 & row > 1) {
                par(mgp = c(0.05, ifelse(row%%2 == 1, 0.15, 0.65),
                  0))
                axis(2, col = axis.col, lwd = 0.5)
            }
            if (col == 1 & row == 1) {
                par(mgp = c(0.05, ifelse(row%%2 == 1, 0.15, 0.65),
                  0))
                axis(2, col = axis.col, lwd = 0.5)
            }
        }
    }
}
<environment: namespace:admbtools>
