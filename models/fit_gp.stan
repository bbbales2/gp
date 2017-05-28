data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  
  vector[2 * N] mu;
  matrix[2 * N, 2 * N] L;
}

transformed data {
  vector[2 * N] z;
  
  z = rep_vector(0, 2 * N);
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
  
  vector[2 * N] yphypph;
}

transformed parameters {
  vector[2 * N] Ay;
  
  for(n in 1:N) {
    Ay[n] = a * y[n] + b * yp[n];
    Ay[n + N] = c * y[n] + d * yp[n];
  }
}

model {
  z ~ multi_normal_cholesky(mu - Ay, L);
  #yphypph ~ multi_normal_cholesky(mu, L);
  
  #yphypph[1:N] ~ normal(a * y + b * yp, sigmay);
  #yphypph[N + 1:2 * N] ~ normal(c * y + d * yp, sigmayp);
  #mu[1:N] ~ normal(a * y + b * yp, sigmay);
  #mu[N + 1:2 * N] ~ normal(c * y + d * yp, sigmayp);
}

