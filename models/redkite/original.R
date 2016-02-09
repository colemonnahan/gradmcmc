marray.juv <- c(42, 18, 5, 7, 4, 3, 2, 1, 2, 2, 1, 0, 1, 3, 0, 0, 1, 1388)
marray.ad <- c(3, 1, 1, 3, 0, 2, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 137)

# Specify model in BUGS language
sink("redkite.jags")
cat("
model {

# Priors and constraints
#sjuv ~ dbeta(4.2, 2.8)   # Informative prior for juv. survival: Analysis A
sjuv ~ dunif(0, 1)      # Non-informative for juv. survival prior: Analysis B
ssub ~ dunif(0, 1)       # Prior for subad. survival
sad ~ dunif(0, 1)        # Prior for ad. survival
rjuv ~ dunif(0, 1)       # Prior for juv. recovery
rad ~ dunif(0, 1)        # Prior for ad. recovery

# Define the multinomial likelihoods
marr.j[1:(n.age+1)] ~ dmulti(pr.j[], rel.j)
marr.a[1:(n.age+1)] ~ dmulti(pr.a[], rel.a)
# Define the cell probabilities of the juvenile m-array
# First element
pr.j[1] <- (1-sjuv)*rjuv
# Second element
pr.j[2] <- sjuv*(1-ssub)*rad
# Third and further elements
for (t in 3:n.age){
   pr.j[t] <- sjuv*ssub*pow(sad,(t-3))*(1-sad)*rad
   }
# Probability of non-recovery
pr.j[n.age+1] <- 1-sum(pr.j[1:n.age])
# Define the cell probabilities of the adult m-array
# All elements
for (t in 1:n.age){
   pr.a[t] <- pow(sad,(t-1))*(1-sad)*rad
   }
# Probability of non-recovery
pr.a[n.age+1] <- 1-sum(pr.a[1:n.age])
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.j = marray.juv, marr.a = marray.ad, n.age = length(marray.juv)-1, rel.j = sum(marray.juv), rel.a = sum(marray.ad))

saveRDS(jags.data, file='data.RDS')

# Initial values
inits <- function(){list(sjuv = runif(1, 0, 1), ssub = runif(1, 0, 1), sad = runif(1, 0, 1), rjuv = runif(1, 0, 1), rad = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("sjuv", "ssub", "sad", "rjuv", "rad")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
rk.ageA <- jags(jags.data, inits, parameters, "redkite.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
rk.ageB <- jags(jags.data, inits, parameters, "redkite.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(rk.ageA, digits = 3)
print(rk.ageB, digits = 3)

par(mfrow = c(2, 3), las = 1)
hist(rk.ageB$BUGSoutput$sims.list$sjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile survival")
hist(rk.ageB$BUGSoutput$sims.list$ssub, breaks = 20, col = "gray", main = "", xlab = "Subadult survival")
hist(rk.ageB$BUGSoutput$sims.list$sad, breaks = 20, col = "gray", main = "", xlab = "Adult survival")
hist(rk.ageB$BUGSoutput$sims.list$rjuv, breaks = 20, col = "gray", main = "", xlab = "Juvenile recovery", xlim = c(0, 0.2))
hist(rk.ageB$BUGSoutput$sims.list$rad, breaks = 20, col = "gray", main = "", xlab = "Adult recovery")

plot(density(rk.ageA$BUGSoutput$sims.list$sjuv), ylim = c(0, 5), lwd = 2, main = "", xlab = "Juvenile survival", las = 1)
points(density(rk.ageB$BUGSoutput$sims.list$sjuv), col = "red", type = "l", lwd = 2)
text(x = 0.5, y = 4.8, "Prior distributions", pos = 4, font = 3)
legend(x = 0.6, y = 4.7, legend = c("U(0,1)", "beta(4.2,2.8)"), lwd = c(2, 2), col = c("black", "red"), bty = "n")

quantile(rk.ageA$BUGSoutput$sims.list$ssub-rk.ageA$BUGSoutput$sims.list$sjuv, prob = c(0.025, 0.975))
quantile(rk.ageA$BUGSoutput$sims.list$sad-rk.ageA$BUGSoutput$sims.list$ssub, prob = c(0.025, 0.975))
