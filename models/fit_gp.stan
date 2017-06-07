data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
}

transformed data {
  vector[N] mu = rep_vector(0, N);
  real ya[N];
  real ypa[N];
  
  for(n in 1 : N) {
    ya[n] = y[n];
    ypa[n] = yp[n];
  }
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmayp;
  real<lower=0> sigmaypp;
  
  vector<lower=0.0>[2] rho;
  real<lower=0.0> alphayp;
  real<lower=0.0> alphaypp;
}

model {
  matrix[N, N] Sigmayp = cov_exp_quad(ya, alphayp, rho[1]) .* cov_exp_quad(ypa, alphayp, rho[2]);
  matrix[N, N] Sigmaypp = cov_exp_quad(ya, alphaypp, rho[1]) .* cov_exp_quad(ypa, alphaypp, rho[2]);
  
  for(n in 1:N) {
    Sigmayp[n, n] = Sigmayp[n, n] + sigmayp + 1e-12;
    Sigmaypp[n, n] = Sigmaypp[n, n] + sigmaypp + 1e-12;
  }
  
  rho ~ gamma(4.0, 4.0);
  alphayp ~ normal(0, 1.0);
  alphaypp ~ normal(0, 1.0);
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);
  
  #ypm ~ multi_normal(mu, Sigmayp);
  #yppm ~ multi_normal(mu, Sigmayp);
  
  yd ~ normal(a * y + b * yp, Sigmayp);
  ypd ~ normal(c * y + d * yp, Sigmaypp);
}

