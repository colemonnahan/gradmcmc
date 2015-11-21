#include <TMB.hpp>
template<class Type>
// bounding function to put uniform priors on parameters
Type boundp(Type x, Type min, Type max){
  return min + (max-min)/(1+exp(-x));
}
vector<Type> boundp2(vector<Type> x, Type min, Type max){
  return min + (max-min)/(1+exp(-x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Ngroup);
  DATA_INTEGER(Nobs);
  DATA_VECTOR(log_yobs);
  DATA_IVECTOR(group);
  DATA_VECTOR(day);
  DATA_VECTOR(spacing);
  PARAMETER(beta2);
  PARAMETER(gamma2);
  PARAMETER_VECTOR(logcpue2);
  PARAMETER(logcpue_mean2);
  PARAMETER(logcpue_sd2);
  PARAMETER_VECTOR(logsigma_obs);
  PARAMETER(sigma_obs_mean2);
  PARAMETER(sigma_obs_sd2);


  // // transform the bounded parameters
  Type logcpue_mean=boundp(logcpue_mean2, Type(-2), Type(5));
  Type logcpue_sd=boundp(logcpue_sd2, Type(.001), Type(5));
  Type sigma_obs_mean=boundp(sigma_obs_mean2, Type(-2), Type(0));
  Type sigma_obs_sd=boundp(sigma_obs_sd2, Type(.001),Type(2));
  Type beta= boundp(beta2, Type(.001),Type(.5));
  Type gamma=boundp(gamma2, Type(.001),Type(.2));
  vector<Type> logcpue=boundp2(logcpue2, Type(-5), Type(5));

  // negative log posterior
  Type nll=0;
  //  nll+=logcpue_mean2+logcpue_sd2+sigma_obs_mean2+sigma_obs_sd2+beta2+gamma2;

  // intermediate variable
  Type ypred;
  // skipping uniform priors for now

  // for(int i=0;i<Ngroup;i++){
  //   nll+=logcpue(i) + logsigma_obs(i);
  // }
  // Loop through random effects and do hyperpriors
  for(int i=0;i<Ngroup;i++){
    nll-=dnorm(logcpue(i), logcpue_mean, logcpue_sd, true);
    nll-=dnorm(logsigma_obs(i), sigma_obs_mean, sigma_obs_sd, true);
  }
  // calculate likelihood of the data
  for(int i=0;i<Nobs;i++){
    ypred= exp(logcpue(group(i)))*exp(-day(i)*gamma)*(1-exp(-beta*spacing(i)));
    nll-=dnorm(log_yobs(i), log(ypred), exp(logsigma_obs(group(i))), true);
  }
  REPORT(logcpue_sd);
  REPORT(sigma_obs_sd);
  REPORT(exp(logsigma_obs(group(0))));
  return nll;
}
