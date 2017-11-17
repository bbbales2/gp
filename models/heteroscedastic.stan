data {
  int<lower=1> N;
  int<lower=1> M;
  real x[N];
  matrix[N, M] y;
}

transformed data {
}

parameters {
  real<lower = 0.0> l;
  real<lower = 0.0> sigmaf;
  vector[N] z1;
  vector[N] z2;
}

transformed parameters {
  vector[N] mu;
  vector[N] sigma;
  
  {
    matrix[N, N] Sigma = cov_exp_quad(x, sigmaf, l);
    matrix[N, N] L;
  
    for(i in 1:N) {
      Sigma[i, i] = Sigma[i, i] + 1e-9;
    }
  
    L = cholesky_decompose(Sigma);
    mu = L * z1;
    sigma = exp(L * z2);
  }
}

model {
  l ~ gamma(4.0, 4.0);
  sigmaf ~ normal(0.0, 1.0);
  z1 ~ normal(0.0, 1.0);
  z2 ~ normal(0.0, 1.0);
  
  for(m in 1:M) {
    y[:, m] ~ normal(mu, sigma);
  }
}
