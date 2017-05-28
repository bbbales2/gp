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
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

model {
  yd ~ normal(a * y[2 : N - 1] + b * yp[2 : N - 1], sigmay);
  ypd ~ normal(c * y[2 : N - 1] + d * yp[2 : N - 1], sigmayp);
}

