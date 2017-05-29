data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
}

transformed data {
  vector[N - 2] yd;
  vector[N - 2] ypd;
  real ya[N];
  real ypa[N];
  
  for(n in 2 : N - 1) {
    yd[n - 1] = (y[n + 1] - y[n - 1]) / (t[n + 1] - t[n - 1]);
    ypd[n - 1] = (yp[n + 1] - yp[n - 1]) / (t[n + 1] - t[n - 1]);
  }
  
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
  
  vector[N - 2] etay;
  vector[N - 2] etayp;
}

transformed parameters {
  vector[N - 2] ypm;
  vector[N - 2] yppm;

  {
    matrix[N - 2, N - 2] Sigma = cov_exp_quad(ya[2:(N - 1)], 1.0, rho[1]) .* cov_exp_quad(ypa[2:(N - 1)], 1.0, rho[2]);
    matrix[N - 2, N - 2] L_Sigma;
    
    for(n in 1:N - 2)
      Sigma[n, n] = Sigma[n, n] + 1e-12;
  
    L_Sigma = cholesky_decompose(Sigma);
    
    ypm = alphayp * L_Sigma * etay;
    yppm = alphaypp * L_Sigma * etayp;
  }
}

model {
  rho ~ gamma(4.0, 4.0);
  alphayp ~ normal(0, 1.0);
  alphaypp ~ normal(0, 1.0);
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);
  
  etay ~ normal(0, 1);
  etayp ~ normal(0, 1);
  
  yd ~ normal(a * y[2 : N - 1] + b * yp[2 : N - 1] + ypm, sigmayp);
  ypd ~ normal(c * y[2 : N - 1] + d * yp[2 : N - 1] + yppm, sigmaypp);
}

