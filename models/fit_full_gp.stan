data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}

parameters {
  real<lower=0.0> l;
  real<lower=0.0> alpha;
  real<lower=0.0> sigma;
  vector[N] zn;
}

transformed parameters {
  vector[N] z;
  
  {
    matrix[N, N] Sigma = cov_exp_quad(x, 1.0, l);
    matrix[N, N] L_Sigma;
    
    for(n in 1:N)
      Sigma[n, n] = Sigma[n, n] + 1e-12;
    
    L_Sigma = cholesky_decompose(Sigma);
    
    z = alpha * L_Sigma * zn;
  }
}

model {
  l ~ gamma(4.0, 4.0);
  alpha ~ normal(0, 1.0);
  
  zn ~ normal(0.0, 1.0);

  y ~ normal(z, sigma);
}

