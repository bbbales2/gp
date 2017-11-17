data {
  int<lower = 1> N;
  real x[N];
  real y[N];
  int<lower = 1> P;
  real lp[P];
}

transformed data {
  matrix[P, P] Sigma_P = cov_exp_quad(lp, 1.0, 1.0);
  matrix[N * N, P] lookup;
  
  {
    matrix[P, N * N] exact;
    for(p in 1:P) {
      matrix[N, N] Sigma = cov_exp_quad(x, 1.0, lp[p]);
      for(n in 1:N) {
        Sigma[n, n] = Sigma[n, n] + 1e-10;
      }
      exact[p, ] = to_row_vector(cholesky_decompose(Sigma));
    }
    
    for(p in 1:P) {
      Sigma_P[p, p] = Sigma_P[p, p] + 1e-10;
    }
    
    lookup = (Sigma_P \ exact)';
  }
}

parameters {
  real<lower = min(lp), upper = max(lp)> l;
  real<lower = 0.0> sigma;
  vector[N] z;
}

transformed parameters {
  vector[N] f;
  
  {
    vector[P] lp_ = to_vector(lp);
    vector[P] Kp = exp(-(l - lp_) .* (l - lp_) / 2.0); //to_vector(cov_exp_quad({ l }, lp, 1.0, 1.0));
    matrix[N, N] L = to_matrix(lookup * Kp, N, N, 1);

    f = L * z;
  }
}

model {
  z ~ normal(0, 1);
  l ~ gamma(4.0, 4.0);
  
  y ~ normal(f, sigma);
}