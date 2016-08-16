 // Taken from https://github.com/stan-dev/example-models/blob/master/BPA/Ch.08/mr_mnl_age3.stan

data {
  int Ndata;
  int Nstage;
  int Nyear;
  int Nplant;
  int year[Ndata];
  int plant[Ndata];
  int stage[Ndata];
  vector[Ndata] Pods;
  int toF[Ndata];
}

parameters {
  real<lower=0.0001,upper=5> yearInterceptSD;
  real<lower=0.0001,upper=5> plantInterceptSD;
  real<lower=0.0001,upper=5> plantSlopeSD;
  vector[Nstage] intercept;
  real slope;
  // Random effect vectors
  vector[Nyear] yearInterceptEffect;
  vector[Nplant] plantInterceptEffect;
  vector[Nplant] plantSlopeEffect;
  }

model {
  vector[Ndata] ypred;
  yearInterceptEffect ~ normal(0, yearInterceptSD);
  plantInterceptEffect ~ normal(0, plantInterceptSD);
  plantSlopeEffect ~ normal(0, plantSlopeSD);
  // Priors
  slope~normal(0, 100);
  intercept~normal(0, 100);
  ypred= intercept[stage] +
    yearInterceptEffect[year] +
    plantInterceptEffect[plant] +
  Pods*slope +
    Pods .* plantSlopeEffect[plant];
  toF ~ bernoulli_logit(ypred);
}

