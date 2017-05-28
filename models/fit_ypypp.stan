data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yph;
  vector[N] ypph;
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

model {
  yph ~ normal(a * y + b * yp, sigmay);
  ypph ~ normal(c * y + d * yp, sigmayp);
}

