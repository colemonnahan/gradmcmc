tits <- read.table("tits.txt", header = TRUE)
C <- as.matrix(tits[5:13])

obs <- as.matrix(tits[14:22])
a <- as.numeric(levels(factor(obs)))     # All the levels, numeric
newobs <- obs                            # Gets ObsID from 1:271
for (j in 1:length(a)) {
  newobs[which(obs == a[j])] <- j
}
newobs[is.na(newobs)] <- 272

first <- as.matrix(tits[23:31])
first[is.na(first)] <- 0

## Separate missing data
df <- expand.grid(site = 1:nrow(C), year = 1:ncol(C))
df$count <- c(C)
obsdata <- subset(df, !is.na(count))
misdata <- subset(df, is.na(count))

## Bundle data. The NAs require different formats for the two platforms
data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), nobs = 272, newobs
             = t(newobs), first = t(first), year = ((1:9)-5) / 4)
data.stan <- list(nobs = nrow(obsdata),
                  nmis = nrow(misdata),
                  nyear = ncol(C),
                  nsite = nrow(C),
                  obs = obsdata$count,
                  obsyear = obsdata$year,
                  obssite = obsdata$site,
                  misyear = misdata$year,
                  missite = misdata$site)

## matplot(1999:2007, t(C), type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Territory counts", xlab = "Year", ylim = c(0, 80), frame = FALSE)
## table(obs)
## length(table(obs))
## apply(first, 2, sum, na.rm = TRUE)

a <- as.numeric(levels(factor(obs)))     # All the levels, numeric
newobs <- obs                            # Gets ObsID from 1:271
for (j in 1:length(a)){newobs[which(obs==a[j])] <- j }
table(newobs)

newobs[is.na(newobs)] <- 272
table(newobs)
first[is.na(first)] <- 0
table(first)


## (h) The full model
## Specify model in BUGS language
sink("tits.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend
beta2 ~ dnorm(0, 0.01)               # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 3)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 1)

for (k in 1:nobs){
   gamma[k] ~ dnorm(0, tau.gamma)   # Random observer effects
   }
tau.gamma <- 1/ (sd.gamma * sd.gamma)
sd.gamma ~ dunif(0, 1)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + gamma[newobs[i,j]] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()


# Initial values
inits.fn <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -1, 1), gamma = runif(272, -1, 1), eps = runif(9, -1, 1))
inits <- list(inits.fn())

# Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps", "sd.alpha", "sd.gamma", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 1

# Call JAGS from R (BRT 11 min)
## out7 <- jags(data, inits, params, "tits.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
## print(out7, dig = 2)

saveRDS(data, 'data.RDS')
saveRDS(data.stan, 'data.stan.RDS')
saveRDS(inits, 'inits.RDS')

