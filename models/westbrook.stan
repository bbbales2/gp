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
  real x[N];
  int y[N];
  real scale;
}

transformed data {
 # real sigma = 1.0;
}

parameters {
  vector[M] z;
  real<lower = 0.0> sigma;
  real<lower = 0.0> l;
}

transformed parameters {
  vector[N] q = approx_L(M, scale, x, sigma, l) * z;
}

model {
  l ~ gamma(4, 4);
  z[1] ~ normal(0.0, 0.01);
  z[2:M] ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ bernoulli_logit(q);
}

generated quantities {
  #vector[N] f;
  real error = -20;
  {
    matrix[N, M] L = approx_L(M, scale, x, sigma, l);
    
    #f = inv_logit(L * z);
    
    error = log10(max(fabs(cov_exp_quad(x, sigma, l) - L * L')) + 1e-20);
  }
}
