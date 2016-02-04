
transformed data {
  real mu;
  real<lower=0> tau;
  real alpha;
  int N;
  mu <- 8;
  tau <- 3;
  N <- 1;
}
parameters {
  real x;
}

model {
  x ~ normal(0, 1);
}

generated quantities {
  real z;
  real t;
  real t2;
  z <- normal_rng(mu, tau);
  t <- student_t_rng(4, mu, tau);
  t2 <- student_t_rng(4, 0,1)*tau+mu;
}

