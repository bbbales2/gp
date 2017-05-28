data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  real alphayp;
  real alphaypp;
  real rho;
}

transformed data {
  vector[N - 2] yd;
  vector[N - 2] ypd;
  
  matrix[N - 2, N - 2] Sigmayp = cov_exp_quad(t[2:(N - 1)], alphayp, rho);
  matrix[N - 2, N - 2] L_Sigmayp;
  matrix[N - 2, N - 2] Sigmaypp = cov_exp_quad(t[2:(N - 1)], alphaypp, rho);
  matrix[N - 2, N - 2] L_Sigmaypp;
  
  for(n in 2 : N - 1) {
    yd[n - 1] = (y[n + 1] - y[n - 1]) / (t[n + 1] - t[n - 1]);
    ypd[n - 1] = (yp[n + 1] - yp[n - 1]) / (t[n + 1] - t[n - 1]);
  }
    
  for(n in 1:N - 2) {
    Sigmayp[n, n] = Sigmayp[n, n] + 1e-8;
    Sigmaypp[n, n] = Sigmaypp[n, n] + 1e-8;
  }
    
  L_Sigmayp = cholesky_decompose(Sigmayp);
  L_Sigmaypp = cholesky_decompose(Sigmaypp);
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmayp;
  real<lower=0> sigmaypp;
  
  vector[N - 2] etay;
  vector[N - 2] etayp;
}

transformed parameters {
  vector[N - 2] ypm;
  vector[N - 2] yppm;
  
  ypm = L_Sigmayp * etay;
  yppm = L_Sigmaypp * etayp;
}

model {
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);
  
  etay ~ normal(0, 1);
  etayp ~ normal(0, 1);
  
  yd ~ normal(a * y[2 : N - 1] + b * yp[2 : N - 1] + ypm, sigmayp);
  ypd ~ normal(c * y[2 : N - 1] + d * yp[2 : N - 1] + yppm, sigmaypp);
}

