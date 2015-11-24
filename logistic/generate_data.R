message('Generating data for logistic model')
setwd('logistic')
## Generate fake data
logistic.traj <- function(logr, logK, num.years, sd.catch, prop.caught,
                          years.obs, sd.obs){
    r <- exp(logr)
    K <- exp(logK)
    catches <- trajectory <- rep(0,num.years)
    trajectory[1] <- K
    for( yr in 2: num.years){
        Nt <- trajectory[yr-1]
        Ct <- catches[yr-1] <-
            if(yr<50) exp(rnorm(1,log(prop.caught*Nt),sd=sd.catch)) else 0
        trajectory[yr] <- Nt+r*Nt*(1- (Nt/K))-Ct
        if( trajectory[yr]<=0 ) { trajectory[yr] <- NA; break}
    }
    ## generate lognormal samples
    pop.obs <-
        exp(rnorm(n=length(years.obs),
                  mean=log(trajectory[years.obs]),sd=sd.obs)
            -sd.obs^2/2)
    ## plot(1:num.years, ylim=c(0,K*1.1),y= trajectory, type='b')
    ## points(1:num.years, catches, type='h')
    ## points(years.obs, pop.obs, col='red')
    return(list(num_years=num.years, catches=catches, num_obs=length(pop.obs), pop_obs=pop.obs,
                sd_obs=sd.obs, years_obs=years.obs))
    ## return(trajectory)
}
logr.true <- log(.1)
logK.true <- log(5000)
set.seed(seeds.list$logistic)
data.logistic <-
    logistic.traj(logr=logr.true, logK=logK.true, num.years=100,
                  sd.catch=.5, prop.caught=.05, years.obs=c(15,35,85,90),
                  sd.obs=.1)
data.logistic.TMB <- data.logistic
saveRDS(data.logistic.TMB, file=results.file("data.logistic.TMB.RData"))
parameters.logistic.TMB <- list(logr=logr.true, logK=logK.true)
saveRDS(parameters.logistic.TMB, file=results.file("parameters.logistic.TMB.RData"))
setwd(main.dir)

message("Generating data for 2d MVN")

