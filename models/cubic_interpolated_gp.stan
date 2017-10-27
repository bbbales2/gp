functions {
  matrix approx_L(real l, real[] lp, matrix[] Ls, matrix[] dLdls);
}


data {
  int<lower = 1> N;
  real x[N];
  real y[N];
  int<lower = 1> P;
  real lp[P];
  matrix[N, N] Ls[P];
  matrix[N, N] dLdls[P];
}

parameters {
  real<lower = min(lp), upper = max(lp)> l;
  real<lower = 0.0> sigma;
  vector[N] z;
}

model {
  vector[N] f;
  f = approx_L(l, lp, Ls, dLdls) * z;
  
  z ~ normal(0, 1);
  l ~ gamma(4.0, 4.0);
  
  y ~ normal(f, sigma);
}