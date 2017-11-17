data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
  real alphayp;
  real alphaypp;
  real rho;
}

transformed data {
  matrix[N, N] Sigma_y;
  matrix[N, N] Sigma_yp;
  matrix[N, N] L_Sigma_y;
  matrix[N, N] L_Sigma_yp;
  
  Sigma_y = cov_exp_quad(to_array_1d(y), 1.0, rho);
  Sigma_yp = cov_exp_quad(to_array_1d(yp), 1.0, rho);
  
  for(n in 1:N) {
    Sigma_y[n, n] = Sigma_y[n, n] + 1e-12;
    Sigma_yp[n, n] = Sigma_yp[n, n] + 1e-12;
  }

  L_Sigma_y = cholesky_decompose(Sigma_y);
  L_Sigma_yp = cholesky_decompose(Sigma_yp);
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmayp;
  real<lower=0> sigmaypp;
  
  vector[N] za;
  vector[N] zb;
  vector[N] zc;
  vector[N] zd;
}

transformed parameters {
  vector[N] aya = alphayp * L_Sigma_y * za;
  vector[N] byb = alphaypp * L_Sigma_yp * zb;
  vector[N] cyc = alphayp * L_Sigma_y * zc;
  vector[N] dyd = alphaypp * L_Sigma_yp * zd;
}

model {
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);

  za ~ normal(0, 1);
  zb ~ normal(0, 1);
  zc ~ normal(0, 1);
  zd ~ normal(0, 1);

  yd ~ normal(a * (y + aya) + b * (yp + byb), sigmayp);
  ypd ~ normal(c * (y + cyc) + d * (yp + dyd), sigmaypp);
}

generated quantities {
  vector[N] lyp;
  vector[N] lypp;
  
  lyp = a * (y + aya) + b * (yp + byb);
  lypp = c * (y + cyc) + d * (yp + dyd);
}
