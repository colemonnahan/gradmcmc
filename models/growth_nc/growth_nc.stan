data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  real loglengths[Nobs];  // observed log cpue
  int fish[Nobs];      // vector to index fish
  int ages[Nobs];      // observed ages
}
parameters {
  // the hyperparameters with uniform priors on them
  real logLinf_mean;
  real<lower=0.01, upper=.5> logLinf_sigma;
  real logk_mean;
  real<lower=0.01, upper=.5> logk_sigma;
  // fixed effects
  real<lower=0, upper=.5> sigma_obs; // data on log scale
  // random effects, standard normal and then trasnformed below
  vector<lower=-5, upper=5>[Nfish] logLinf_raw;
  vector<lower=-5, upper=5>[Nfish] logk_raw;
  real<lower=0, upper=5> delta;
}

transformed parameters{
 vector[Nfish] logLinf;
 vector[Nfish] logk;
 // Matt trick: implies logLinf~N(logLinf_mean, logLinf_sigma); etc.
 logLinf  <- logLinf_mean+logLinf_sigma*logLinf_raw;
 logk <- logk_mean+logk_sigma*logk_raw;
 }
model {
 // Loop through random effects and do hyperpriors
 real Linf;
 real k;
 real ypred[Nobs];
 // hyperparams
 logLinf_raw~normal(0,1);
 logk_raw~normal(0,1);

 // calculate likelihood of data
 for(i in 1:Nobs){
  Linf <- exp(logLinf[fish[i]]);
  k <- exp(logk[fish[i]]);
  ypred[i] <- log( Linf*(1-exp(-k*(ages[i]-5)))^delta );
 }
  loglengths~normal(ypred, sigma_obs);
}