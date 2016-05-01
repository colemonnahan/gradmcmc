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
  real sa2;
  real sm2;
  real st2;
  real sf2;
  real se2;
  // random effects
  real ua[n_animal];
  real um[n_ID];
  real ut[n_Tank_ID];
  real uf[n_Dam];
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
  ua~normal(0, sqrt(sa2));
  um~normal(0, sqrt(sm2));
  ut~normal(0, sqrt(st2));
  uf~normal(0, sqrt(sf2));

// model predictions and likelihood
for(i in 1:N){
  ypred[i] <- a[cross[i]]+ b*days[i]+ ua[animal[i]]+ um[ID[i]]+
    ut[Tank_ID[i]] +uf[Dam[i]];
}
  // likelihood
  body_weight~normal(ypred, sqrt(se2));

}


