data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
}

transformed data {
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
  
  vector[N] etay;
  vector[N] etayp;
}

transformed parameters {
  vector[N] ypm;
  vector[N] yppm;

  {
    matrix[N, N] Sigma = cov_exp_quad(ya, 1.0, rho[1]) .* cov_exp_quad(ypa, 1.0, rho[2]);
    matrix[N, N] L_Sigma;
    
    for(n in 1:N)
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
  
  yd ~ normal(a * y + b * yp + ypm, sigmayp);
  ypd ~ normal(c * y + d * yp + yppm, sigmaypp);
}

