// Simple logistic model

data {
  int<lower=0> num_years;
  real catches[num_years];
  int<lower=0> num_obs;
  real log_pop_obs[num_obs];
  real sd_obs;
  int years_obs[num_obs];
  }
parameters {
  real<lower=log(.01), upper=log(.5)> logr;
  real<lower=log(1000), upper=log(15000)> logK;
}

transformed parameters {
 real r;
 real K;
 r <- exp(logr);
 K <- exp(logK);
}
model {
  real N[num_years];
  real ypred[num_obs];
  real temp;
  N[1] <- K;
  ## Loop through and project population dynamics
  for(i in 2:num_years){
        temp <- N[i-1]+r*N[i-1]*(1-N[i-1]/K)-catches[i];
     if(temp<1){
	// increment_log_prob(-.0001*(temp-10)^2);
	N[i] <- 10/(2-temp/10);
	// if(i==20) print(temp,  N[i],  -.0001*(temp-10)^2);
	} else {
	N[i] <- temp;
     }
  }
 // calculate likelihood of data
 for(i in 1:num_obs){
  ypred[i] <- log(N[years_obs[i]]);
 }
  log_pop_obs~normal(ypred, sd_obs);
}