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
  matrix[N, N] K = cov_exp_quad(t, alpha, rho);
  matrix[N, N] L_K;
  real sq_sigma = square(sigma);
  
  // diagonal elements
  for (n in 1:N) K[n, n] = K[n, n] + sq_sigma;
  
  L_K = cholesky_decompose(K);
  
  rho ~ gamma(4, 4);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y ~ multi_normal_cholesky(mu, L_K);
}
