// The Hamley and Skud effective hook model.

data {
  int<lower=0> Npar;
  matrix[Npar,Npar] covar;
  vector[Npar] x;
}
parameters {
 vector[Npar] mu;
}

model {
  x~multi_normal(mu, covar);
}