data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
}

transformed data {
  vector[N - 2] yd;
  vector[N - 2] ypd;
  
  for(n in 2 : N - 1) {
    yd[n - 1] = (y[n + 1] - y[n - 1]) / (t[n + 1] - t[n - 1]);
    ypd[n - 1] = (yp[n + 1] - yp[n - 1]) / (t[n + 1] - t[n - 1]);
  }
}

parameters {
  real b;
  real c;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

model {
  yd ~ normal(b * yp[2 : N - 1], sigmay);
  ypd ~ normal(c * sin(y[2 : N - 1]), sigmayp);
}

