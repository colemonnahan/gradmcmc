// The Hamley and Skud effective hook model.

data {
  int<lower=0> Ngroup; // number of groups
  int<lower=0> Nobs;   // number of observations
  real log_yobs[Nobs];  // observed log cpue
  int group[Nobs];      // vector to index group
  int day[Nobs];       // day effect (should this be integer or real?)
  real spacing[Nobs];   // hook spacing (indepedent variable)
}
parameters {
  // the hyperparameters with uniform priors on them
  real<lower=-5, upper=5> logcpue_mean;
  real<lower=0, upper=5> logcpue_sd;
  real<lower=-1, upper=5> sigma_obs_mean;
  real<lower=0, upper=3> sigma_obs_sd;
  // fixed effects
  real<lower=0, upper=10> beta;
  real<lower=0, upper=1> gamma;
  // random effects
  real logcpue[Ngroup];
  real logsigma_obs[Ngroup];
}
// transformed parameters {
//   for (j in 1:J)
//     theta[j] <- mu + tau * eta[j];
// }
model {

// explicit priors even though they are uniform
logcpue_mean~uniform(-5,5);
logcpue_sd~uniform(0,5);
sigma_obs_mean~uniform(-1,5);
sigma_obs_sd~uniform(0,3);
beta~uniform(0,10);
gamma~uniform(0,1);

// Loop through random effects and do hyperpriors
 for(i in 1:Ngroup){
 // mean cpue for each site in log space
 logcpue[i]~normal(logcpue_mean, logcpue_sd);
 logsigma_obs[i]~normal(sigma_obs_mean, sigma_obs_sd);
 }

 // calculate likelihood of data
 for(i in 1:Nobs){
 real ypred;
  ypred <- exp(logcpue[group[i]])*exp(-day[i]*gamma)*(1-exp(-beta*spacing[i]));
  //ypred <- logcpue[group[i]]+(-day[i]*gamma)+(1-exp(-beta*spacing[i]));
 //if(i==1) {print(ypred); print(beta);}
 log_yobs[i]~normal(log(ypred), exp(logsigma_obs[group[i]]));
 }
}