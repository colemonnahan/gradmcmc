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
  real<lower=0> yearInterceptSD;
  real<lower=0> plantInterceptSD;
  real<lower=0> plantSlopeSD;
  vector[Nstage] intercept;
  real slope;
  // Random effect vectors
  vector[Nyear] yearInterceptEffect;
  vector[Nplant] plantInterceptEffect;
  vector[Nplant] plantSlopeEffect;
}

model {
  vector[Ndata] ypred;
  // Random effects, centered
  yearInterceptEffect ~ normal(0, yearInterceptSD);
  plantInterceptEffect ~ normal(0, plantInterceptSD);
  plantSlopeEffect ~ normal(0, plantSlopeSD);
  // Priors
  yearInterceptSD~cauchy(0,5);
  plantInterceptSD~cauchy(0,5);
  plantSlopeSD~cauchy(0,5);
  slope~normal(0, 10);
  intercept~normal(0, 10);
  // vectorized prediction
  ypred= intercept[stage] +
    yearInterceptEffect[year] +
    plantInterceptEffect[plant] +
    Pods*slope + Pods .* plantSlopeEffect[plant];
  toF ~ bernoulli_logit(ypred);
}

