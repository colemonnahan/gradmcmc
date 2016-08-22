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
  vector[Nyear] yearInterceptEffect_raw;
  vector[Nplant] plantInterceptEffect_raw;
  vector[Nplant] plantSlopeEffect_raw;
}

model {
  vector[Ndata] ypred;
  // noncentered random effects; implies N(0,SD) below
  yearInterceptEffect_raw ~ normal(0, 1);
  plantInterceptEffect_raw ~ normal(0, 1);
  plantSlopeEffect_raw ~ normal(0, 1);
  // Priors
  yearInterceptSD~cauchy(0,5);
  plantInterceptSD~cauchy(0,5);
  plantSlopeSD~cauchy(0,5);
  slope~normal(0, 10);
  intercept~normal(0, 10);
  ypred= intercept[stage] +
    yearInterceptEffect_raw[year]*yearInterceptSD +
    plantInterceptEffect_raw[plant]*plantInterceptSD +
    Pods .* (plantSlopeEffect_raw[plant]*plantSlopeSD)+
    Pods*slope;
  toF ~ bernoulli_logit(ypred);
}

