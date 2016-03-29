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
  // random effects, needed to boudn these since they were going crazy during tuning
  real<lower=2, upper=5> logLinf[Nfish];
  real<lower=-5,upper=-1> logk[Nfish];
  real<lower=0, upper=5> delta;
}

model {

 // Loop through random effects and do hyperpriors
 real Linf;
 real k;
 real ypred[Nobs];
 // hyperparams
 logLinf~student_t(100, logLinf_mean, logLinf_sigma);
 logk~student_t(100,logk_mean, logk_sigma);

 // calculate likelihood of data
 for(i in 1:Nobs){
  Linf <- exp(logLinf[fish[i]]);
  k <- exp(logk[fish[i]]);
  ypred[i] <- log( Linf*(1-exp(-k*(ages[i]-5)))^delta );
 }
  loglengths~student_t(100, ypred, sigma_obs);
}