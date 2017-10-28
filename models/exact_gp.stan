data {
  int<lower = 1> N;
  real x[N];
  real y[N];
}

parameters {
  real<lower = 0.0> l;
  real<lower = 0.0> sigma;
  vector[N] z;
}

transformed parameters {
  vector[N] f;
  
  {
    matrix[N, N] Sigma = cov_exp_quad(x, 1.0, l);
    matrix[N, N] L;
    
    for(n in 1:N) {
      Sigma[n, n] = Sigma[n, n] + 1e-10;
    }
    L = cholesky_decompose(Sigma);
    
    f = L * z;
  }
}

model {
  z ~ normal(0, 1);
  l ~ gamma(4.0, 4.0);
  
  y ~ normal(f, sigma);
}