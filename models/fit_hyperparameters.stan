data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
}

transformed data {
  vector[N] mu;
  mu = rep_vector(0, N);
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N, N] Sigma = cov_exp_quad(t, alpha, rho);
  matrix[N, N] L_Sigma;
  real sq_sigma = square(sigma);
  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = Sigma[k, k] + sq_sigma;
  L_Sigma = cholesky_decompose(Sigma);
  
  rho ~ gamma(4, 4);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);
}

