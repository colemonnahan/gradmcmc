/*
 * A type of Cormack-Jolly-Seber (CJS) model
 *
 *
 * References
 *
 *    Cormack, R. M. 1964. Estimates of survival from the sighting of
 *         marked animals. Biometrika 51:429-438.)
 *    Seber, G. A. F. 1965. A note on the multiple-recapture
 *         census. Biometrika 52:249-259.
 *
 * see also:
 *
 *    Lebreton, J.-D., K. P. Burnham, J. Clobert, and
 *    D. R. Anderson. 1992.  Modelling survival and testing biological
 *    hypotheses using marked animals: a unified approach with case
 *    studies. Ecological Monographs 62:67-118.
 *
 * for BUGS code, see:
 *
 *    Kéry, M., and M. Schaub. 2012. Bayesian Population Analysis
 *         using WinBUGS. Elsevier, Amsterdam.
 *
 * Specific characteristics of the CJS-model coded below:
 *
 *  - all I individuals are released at time k=0
 *  - capture probability p depends on a linear predictor
 *    including a random family effect
 *  - survival probability phi depends on a linear predictor
 *    including a random family effect
 *
 */

data {
  int<lower=2> K;                       // capture events
  int<lower=0> I;                       // number of individuals
  int<lower=0,upper=1> CH[I,K];         // CH[i,k]: individual i captured at k
  int<lower=0> nfam;                    // number of families
  int<lower=0, upper=nfam> family[I];   // index of group variable
  vector[I] carez;                     // covariable, duration of parental care, z-trans.
  int<lower=1,upper=4> year[I];         // index of year
  vector[K] agec;                      // age of fledling, centered
}

transformed data {
  int<lower=0,upper=K+1> last[I];       // last[i]:  ind i last capture
  last <- rep_array(0,I);
  for (i in 1:I) {
    for (k in 1:K) {
      if (CH[i,k] == 1) {
        if (k > last[i])  last[i] <- k;
      }
    }
  }
}

parameters {
  real b0[4];                  // intercepts per year for p
  real b1[4];                  // slope for age per year for p
  real a[K-1];                   // intercept of phi
  real a1;                     // coef of phi
  real<lower=0> sigmaphi;      // between family standard deviation in logit(phi)
  real<lower=0> sigmayearphi;  // between-year standard deviation in logit(phi)
  real<lower=0> sigmap;        // between family standard deviation in logit(p)
  real fameffphi[nfam];        // family effects for phi
  real fameffp[nfam];          // family effects for p
  real yeareffphi[4];          // year effect on phi
}

transformed parameters {
  real<lower=0,upper=1>p[I,K];       // capture probability
  real<lower=0,upper=1>phi[I,K-1];   // survival probability
  real<lower=0,upper=1>chi[I,K+1];   // probability that an individual
                                     // is never recaptured after its
                                     // last capture
  {
    int k;
    for(ii in 1:I){
      for(tt in 1:(K-1)) {
        // linear predictor with random effect for phi:
        // add fixed and random effects here
        phi[ii,tt] <- inv_logit(a[tt] +a1*carez[ii] +
          sigmayearphi * yeareffphi[year[ii]] + sigmaphi*fameffphi[family[ii]]);
      }
    }
    for(i in 1:I) {
      // linear predictor with random effect
      // for p: add fixed and random effects here
      p[i,1] <- 1;  // first occasion is marking occasion
      for(kk in 2:K)
        p[i,kk] <- inv_logit(b0[year[i]] + b1[year[i]]*agec[kk]+
                                sigmap*fameffp[family[i]]);
      chi[i,K+1] <- 1.0;
      k <- K;
      while (k > 1) {
        chi[i,k] <- (1 - phi[i,k-1]) + phi[i,k-1] * (1 - p[i,k]) * chi[i,k+1];
        k <- k - 1;
      }
      chi[i,1] <- (1 - p[i,1]) * chi[i,2];
    }
  }
}


model {
  // priors
  for(j in 1:4){
      b0[j]~normal(0, 10);
      b1[j]~normal(0, 10);
      }
  for(v in 1:(K-1)){
    a[v]~normal(0,10);
  }
  a1~normal(0,10);
  sigmaphi~uniform(0,5);
  sigmayearphi~uniform(0,3);
  sigmap~uniform(0,5);

  // random effects
  for(g in 1:nfam) {
    fameffphi[g]~normal(0, 1);
    fameffp[g]~normal(0,1);
  }
  for(ye in 1:4){
      yeareffphi[ye]~normal(0, 1);
  }
  // likelihood
  for (i in 1:I) {
    if (last[i]>0) {
      for (k in 1:last[i]) {
        if(k>1) increment_log_prob(log(phi[i, k-1]));
        if (CH[i,k] == 1) increment_log_prob(log(p[i,k]));
        else increment_log_prob(log1m(p[i,k]));
      }
    }
    increment_log_prob(log(chi[i,last[i]+1]));
  }
}
