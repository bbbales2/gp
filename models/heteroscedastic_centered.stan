data {
  int<lower=1> N;
  int<lower=1> M;
  real x[N];
  matrix[N, M] y;
}

transformed data {
  vector[N] zeros = rep_vector(0.0, N);
}

parameters {
  real<lower = 0.0> l;
  real<lower = 0.0> sigmaf;
  vector[N] mu;
  vector<lower = 0.0>[N] sigma_log;
}

transformed parameters {
  vector[N] sigma = exp(sigma_log);
}

model {
  matrix[N, N] Sigma = cov_exp_quad(x, sigmaf, l);
  matrix[N, N] L;
  
  for(i in 1:N) {
    Sigma[i, i] = Sigma[i, i] + 1e-9;
  }
  
  L = cholesky_decompose(Sigma);
  
  mu ~ multi_normal_cholesky(zeros, L);
  sigma_log ~ multi_normal_cholesky(zeros, L);
  l ~ gamma(4.0, 4.0);
  sigmaf ~ normal(0.0, 1.0);
  
  for(m in 1:M) {
    y[:, m] ~ normal(mu, sigma);
  }
}
