// The Hamley and Skud effective hook model.

data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  real lengths[Nobs];  // observed log cpue
  int fish[Nobs];      // vector to index fish
  int ages[Nobs];      // observed ages
}
parameters {
  // the hyperparameters with uniform priors on them
  real<lower=30, upper=70> Linf_mean;
  real<lower=0, upper=.05> Linf_sigma;
  real<lower=0, upper=.5> k_mean;
  real<lower=0, upper=.05> k_sigma;
  // fixed effects
  real<lower=0, upper=1> sigma_obs; // data on log scale
  // random effects
  real<lower=10, upper=100> Linf[Nfish];
  real<lower=0,upper=.5> k[Nfish];
}
// transformed parameters {
//   for (j in 1:J)
//     theta[j] <- mu + tau * eta[j];
// }
model {
// explicit priors even though they are uniform
 // ??

// Loop through random effects and do hyperpriors
 for(i in 1:Nfish){
 // mean cpue for each site in log space
 Linf[i]~normal(Linf_mean, Linf_sigma);
 k[i]~normal(k_mean, k_sigma);
 }

 // calculate likelihood of data
 for(i in 1:Nobs){
 real ypred;
  ypred <- Linf[fish[i]]*(1-exp(-k[fish[i]]*(ages[i])));
 //if(i==1) {print(ypred); print(beta);}
 lengths[i]~normal(log(ypred), sigma_obs);
 }
}