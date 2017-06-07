data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
}

parameters {
  real c;
  real<lower=0> sigmayp;
}

model {
  ypd ~ normal(c * sin(y), sigmayp);
}

