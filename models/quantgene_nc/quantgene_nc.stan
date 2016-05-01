## This is the glmmMCMC code:
## body_weight ~ Reciprocal_Cross + degree_days_centered, random =
##              ~animal + Dam + Tank_ID + ID


data {
  int<lower=0> N;		// # of observations
  vector[N] body_weight;	//
  vector[N] days;
  int<lower=0> cross[N];
  int<lower=0> animal[N];
  int<lower=0> Dam[N];
  int<lower=0> Tank_ID[N];
  int<lower=0> ID[N];
  int n_cross;
  int n_animal;
  int n_Dam;
  int n_Tank_ID;
  int n_ID;
}

parameters {
  // fixed effects
  real a[n_cross];
  real b;
  // hyper variances
  real<lower=0, upper=10> sa2;
  real<lower=0, upper=10> sm2;
  real<lower=0, upper=10> st2;
  real<lower=0, upper=10> sf2;
  real<lower=0, upper=10> se2;
  // random effects
  real ua_raw[n_animal];
  real um_raw[n_ID];
  real ut_raw[n_Tank_ID];
  real uf_raw[n_Dam];
}

transformed parameters {
  // non-centered random effects; implies ua~N(0, sa); etc.
  // real ua[n_animal];
  // real um[n_ID];
  // real ut[n_Tank_ID];
  // real uf[n_Dam];
  // vector[n_animal] ua;
  // vector[n_ID] um;
  // vector[n_Tank_ID] ut;
  // vector[n_Dam] uf;
  // ua <- ua_raw*sqrt(sa2);
  // um <- um_raw*sqrt(sm2);
  // ut <- ut_raw*sqrt(st2);
  // uf <- uf_raw*sqrt(sf2);
  real sa;
  real sm;
  real st;
  real sf;
  sa <- sqrt(sa2);
  sm <- sqrt(sm2);
  st <- sqrt(st2);
  sf <- sqrt(sf2);
}

model {
  real ypred[N];

  // Priors:
  // Fixed effects with vague normal priors
  a~normal(0, sqrt(1e10));		// cross effects
  b~normal(0, sqrt(1e10));		// days effect
  // Random effects: variance terms; inverse Wishart priors
  sa2~student_t(2,0,1);		// additive genetic effects (animal?)
  sm2~student_t(2,0,1);		// maternal effects (ID?)
  st2~student_t(2,0,1);		// tank effects (Tank_ID?)
  sf2~student_t(2,0,1);		// family effects (Dam?)
  se2~student_t(2,0,1);		// residual error

  // Hyperdistribution
  ua_raw~normal(0,1);
  um_raw~normal(0,1);
  ut_raw~normal(0,1);
  uf_raw~normal(0,1);

  // model predictions and likelihood
  for(i in 1:N){
    ypred[i] <- a[cross[i]]+ b*days[i]+ ua_raw[animal[i]]*sa+ um_raw[ID[i]]*sm+ ut_raw[Tank_ID[i]]*st +uf_raw[Dam[i]]*sf;
  }
  // likelihood
  body_weight~normal(ypred, sqrt(se2));

}


