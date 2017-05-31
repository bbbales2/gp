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
  real ya[N];
  real ypa[N];
  matrix[N, N] Sigma;
  matrix[N, N] L_Sigma;
  
  for(n in 1 : N) {
    ya[n] = y[n];
    ypa[n] = yp[n];
  }
  
  Sigma = cov_exp_quad(ya, 1.0, rho) .* cov_exp_quad(ypa, 1.0, rho);
  
  for(n in 1:N)
    Sigma[n, n] = Sigma[n, n] + 1e-12;

  L_Sigma = cholesky_decompose(Sigma);
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmayp;
  real<lower=0> sigmaypp;
  
  vector[N] etay;
  vector[N] etayp;
}

transformed parameters {
  vector[N] ypm;
  vector[N] yppm;
  
  ypm = alphayp * L_Sigma * etay;
  yppm = alphaypp * L_Sigma * etayp;
}

model {
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);
  
  etay ~ normal(0, 1);
  etayp ~ normal(0, 1);
  
  yd ~ normal(a * y + b * yp + ypm, sigmayp);
  ypd ~ normal(c * y + d * yp + yppm, sigmaypp);
}

generated quantities {
  vector[N] lyp;
  vector[N] lypp;
  
  lyp = a * y + b * yp;
  lypp = c * y + d * yp;
}
