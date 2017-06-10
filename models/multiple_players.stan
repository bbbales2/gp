functions {
  matrix approx_L(int M, real scale, real[] xt, real sigma, real l) {
    int N = size(xt);
    real a = 1.0 / (4.0 * scale^2);
    real b = 1.0 / (2.0 * l^2);
    real c = sqrt(a^2.0 + 2.0 * a * b);
    
    real epsilon = sqrt(b);
    real alpha = sqrt(2.0 * a);
    real beta = (1.0 + (2.0 * epsilon / alpha)^2)^0.25;
    real delta = sqrt(alpha^2 * (beta^2 - 1.0) / 2.0);
    
    vector[N] Ht[M];
    vector[N] x = to_vector(xt);
    matrix[N, M] Htt;
    vector[N] xp = alpha * beta * x;
    real f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2));
    
    Ht[1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x .* x);
    Ht[2] = f * sqrt(1.0 / 2) * 2.0 * xp .* Ht[1];
    for(n in 3:M) {
      Ht[n] = f * sqrt(1.0 / (2.0 * (n - 1))) * 2.0 * xp .* Ht[n - 1] - f^2 * sqrt(1.0 / (4.0 * (n - 1) * (n - 2))) * 2.0 * (n - 2) * Ht[n - 2];
    }
    
    for(n in 1:M) {
      Htt[:, n] = Ht[n];
    }
    
    return sigma * Htt;
  }
}

data {
  int N;
  int M;
  int L;
  real x[N * L];
  int y[N * L];
  real scale;
}

parameters {
  vector[M] zt;
  real<lower = 0.0> sigmat;
  real<lower = 0.0> lt;
  vector[M] z[L];
  real<lower = 0.0> sigma[L];
  real<lower = 0.0> l[L];
}

transformed parameters {
  vector[N] f[L];
  vector[N * L] mu = approx_L(M, scale, x, sigmat, lt) * zt;
  for(ll in 1:L) {
    f[ll] = inv_logit(mu[(ll - 1) * N + 1 : ll * N] + approx_L(M, scale, x[(ll - 1) * N + 1 : ll * N], sigma[ll], l[ll]) * z[ll]);
  }
}

model {
  lt ~ gamma(4, 4);
  zt ~ normal(0, 1);
  sigmat ~ normal(0, 1);
  for(ll in 1:L) {
    z[ll] ~ normal(0, 1);
    sigma[ll] ~ normal(0, 1);
    l[ll] ~ gamma(4, 4);
  }
  
  for(ll in 1:L) {
    y[(ll - 1) * N + 1 : ll * N] ~ bernoulli(f[ll]);
  }
}