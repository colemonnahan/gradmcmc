#include <TMB.hpp>
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(num_years);
  DATA_VECTOR(catches);
  DATA_INTEGER(num_obs);
  DATA_VECTOR(pop_obs);
  DATA_SCALAR(sd_obs);
  DATA_VECTOR(years_obs);
  PARAMETER(logr);
  PARAMETER(logK);
  vector<Type> N(num_years);
  int j=0;			// counter
  N(0)=exp(logK);
  Type pen=0; // penalty for crashing population (posfun)
  Type eps=1;
  Type nll=0; 			// negative log likelihood
  vector<Type> temp(num_years);
  for(int t=1; t<num_years; t++){
    temp(t)=N(t-1)+exp(logr)*N(t-1)*(1-N(t-1)/exp(logK))-catches(t-1);
    // use posfun to keep from going negative
    N(t)=posfun(temp(t), eps, pen);
    if(j<num_obs){
      if(years_obs(j)==t){
	nll-=dnorm(log(pop_obs(j)), log(N(t)), sd_obs, true);
	j++;
      }
    }
  }
  // add priors
  // nll-=dnorm(exp(logK), Type(5000), Type(1000));
  // nll-=dnorm(exp(logr), Type(.1), Type(.02));
  REPORT(N);
  REPORT(pen);
  REPORT(nll);
  return nll+pen;
}

