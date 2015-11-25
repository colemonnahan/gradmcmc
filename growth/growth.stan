// The Hamley and Skud effective hook model.

data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  real loglengths[Nobs];  // observed log cpue
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
  real<lower=2, upper=5> logLinf[Nfish];
  real<lower=-5,upper=-1> logk[Nfish];
}
model {
// Loop through random effects and do hyperpriors
 for(i in 1:Nfish){
 // mean cpue for each site in log space
 logLinf[i]~normal(log(Linf_mean), Linf_sigma);
 logk[i]~normal(log(k_mean), k_sigma);
 }

 // calculate likelihood of data
 for(i in 1:Nobs){
 real ypred;
  ypred <- logLinf[fish[i]] + log(1-exp(-exp(logk[fish[i]])*ages[i]));
 //if(i==1) {print(ypred); print(beta);}
 loglengths[i]~normal(ypred, sigma_obs);
 }
}