data {
  int<lower=1> N;
  real t[N];
  vector[N] yd;
  vector[N] ypd;
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
  yd ~ normal(a * y + b * yp, sigmay);
  ypd ~ normal(c * y + d * yp, sigmayp);
}

