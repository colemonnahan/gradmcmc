data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  real loglengths[Nobs];  // observed log cpue
  int fish[Nobs];      // vector to index fish
  int ages[Nobs];      // observed ages
}
parameters {
  // fixed effects
  real<lower=0, upper=5> delta;
  real<lower=0, upper=2> sigma_obs; // data on log scale

  // hyperparameters with bounds
  real<lower=-5, upper=5> logLinf_mean;
  real<lower=-5, upper=5> logk_mean;
  real<lower=0, upper=2> logLinf_sigma;
  real<lower=0, upper=2> logk_sigma;

  // random effects
  real<lower=-10, upper=10> logLinf[Nfish];
  real<lower=-10, upper=10> logk[Nfish];
}

model {
  real Linf;
  real k;
  real ypred[Nobs];

  // priors
  // delta is uniform above
  sigma_obs~cauchy(0,2);

  // hyperpriors
  logLinf_sigma~cauchy(0,2);
  logk_sigma~cauchy(0,2);
  // hyper means are uniform above

  // random effects
  logLinf~normal(logLinf_mean, logLinf_sigma);
  logk~normal(logk_mean, logk_sigma);

  // calculate likelihood of data
   for(i in 1:Nobs){
    Linf <- exp(logLinf[fish[i]]);
    k <- exp(logk[fish[i]]);
    ypred[i] <- log( Linf*(1-exp(-k*(ages[i]-5)))^delta );
   }
  loglengths~normal(ypred, sigma_obs);
}
