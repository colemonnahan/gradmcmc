// JAGS model
/* model <- function(){  ## The 8-schools problem in JAGS
    ## priors for fixed effects
    mu~dunif(-20,20)
    tau~dunif(0,100)
    ## Loop through observations and calculate likelihood
    for(i in 1:8){
        ## Hyperdistribution for eta
        eta[i]~dnorm(0,1)
        ## The predicted values
        theta[i] <- mu+ tau*eta[i]
        ## Likelihood of data
        y[i]~dnorm(theta[i], pow(sigma[i],-2))
    }
}*/
#include <TMB.hpp>
template<class Type>
Type boundp(Type x, Type min, Type max){
  return min + (max-min)/(1+exp(-x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(sigma);
  PARAMETER(mu);
  PARAMETER(tau_unbounded);
  PARAMETER_VECTOR(eta);
  Type tau=boundp(tau_unbounded, Type(0), Type(50));
  Type nll=0;
  vector<Type> theta(8);
  for(int i=0;i<=7;i++){
    theta(i)=mu+tau*eta(i);
  // likelihood
    nll-=dnorm(Y(i),theta(i),sigma(i),true);
  // hyperdistribution
    nll-=dnorm(eta(i), Type(0), Type(1), true);
  }
  // add uniform priors to get the log density to match JAGS
  nll+=log(40); // for mu~Unif(-20,20)
  nll+=log(100); // for tau~Unif(0,100)
  REPORT(nll);
  REPORT(theta);
  return nll;
}
// Stan model
/* data {
  int<lower=0> J; // number of schools
  real y[J]; // estimated treatment effects
  real<lower=0> sigma[J]; // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] <- mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
 */
