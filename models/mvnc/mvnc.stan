data {
  int<lower=0> Npar;
  matrix[Npar,Npar] covar;
  vector[Npar] x;
}
parameters {
  vector<lower=-10, upper=10>[Npar] mu;
}

model {
  x~multi_normal(mu, covar);
}
