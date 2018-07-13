data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
}

parameters {
  //real a;
  //real b;
  real c;
  //real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

model {
  yd ~ normal(0 * y + 1.0 * yp, sigmay);
  ypd ~ normal(c * y + 0 * yp, sigmayp);
}

