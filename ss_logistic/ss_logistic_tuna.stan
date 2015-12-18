// The Hamley and Skud effective hook model.

data {
  int<lower=0> N; // number of years
  real catches[N];
  real logcpue[N];
}
parameters {

  real<lower=10, upper=1000> logK;
  real<lower=log(.01), upper=log(1.2)> logr;
  real<lower=.5, upper=100> iq;
  real isigma2;
  real itau2;
  real u[N];
}

transformed parameters {

  real<lower=0> sigma2;
  sigma2 <- 1/isigma2;
  real<lower=0> tau2;
  tau2 <- 1/itau2;
  real<lower=0> q;
  q <- 1/iq;
  real K;
  real r;
  K <- exp(logK);
  r <- exp(logr);
}

model {
 real B[N];
 u[1]~dnorm(0, sqrt(sigma2));
 B[1] <- K
 for(i in 1:N){
 u[i]~dnorm(0, sqrt(sigma2));
 B[i] <- (B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catches[i-1])*exp(u[i])
 ypred[i] <- log(B[i]) +log(*q);
 }
  logcpue~normal(ypred, sqrt(tau2));
}