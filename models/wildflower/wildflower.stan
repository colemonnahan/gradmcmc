 // Taken from https://github.com/stan-dev/example-models/blob/master/BPA/Ch.08/mr_mnl_age3.stan

data {
  int Ndata;
  int Nstage;
  int Nyear;
  int Nplant;
  int year[Ndata];
  int plant[Ndata];
  int stage[Ndata];
  int Pods[Ndata];
  int toF[Ndata];
}

parameters {
  real<lower=0.0001,upper=5> yearInterceptSD;
  real<lower=0.0001,upper=5> plantInterceptSD;
  real<lower=0.0001,upper=5> plantSlopeSD;
  real intercept[Nstage];
  real slope;
  // Random effect vectors
  real yearInterceptEffect[Nyear];
  real plantInterceptEffect[Nplant];
  real plantSlopeEffect[Nplant];
  }

model {
  yearInterceptEffect ~ normal(0, yearInterceptSD);
  plantInterceptEffect ~ normal(0, plantInterceptSD);
  plantSlopeEffect ~ normal(0, plantSlopeSD);
  // Priors
  slope~normal(0, 100);
  intercept~normal(0, 100);
  for(i in 1:Ndata) {
  toF[i]~ bernoulli_logit(
   intercept[stage[i]] + yearInterceptEffect[year[i]] +
   plantInterceptEffect[plant[i]] + slope*Pods[i] +
   plantSlopeEffect[plant[i]] * Pods[i]);
 }
}

