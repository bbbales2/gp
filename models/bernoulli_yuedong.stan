functions {
  real k0(real x) {
    return 1.0;
  }

  real k1(real x) {
    return x - 0.5;
  }
  
  real k2(real x) {
    real k1_ = k1(x);
    return 0.5 * (k1_ * k1_ - 1.0 / 12.0);
  }
  
  real k4(real x) {
    real k1_ = k1(x);
    return (1.0 / 24.0) * (k1_ * k1_ * k1_  * k1_ - 0.5 * k1_ * k1_ + 7.0 / 240.0);
  }
}

data {
  int N;
  real x[N];
  vector[N] y;
}

transformed data {
  cov_matrix[N] Sigma;
  matrix[N, N] L;
  vector[N] vk0;
  vector[N] vk1;
  
  for(i in 1:N) {
    for(j in 1:N) {
      Sigma[i, j] = k2(x[i]) * k2(x[j]) - k4(fabs(x[i] - x[j]));
    }
  }
  
  for(i in 1:N) {
    Sigma[i, i] = Sigma[i, i] + 1e-12;
    
    vk0[i] = k0(x[i]);
    vk1[i] = k1(x[i]);
  }
  
  L = cholesky_decompose(Sigma);
}

parameters {
  real a;
  real b;
  real<lower=0.0> alpha;
  real<lower=0.0> sigma;
  
  vector[N] z;
}

transformed parameters {
  vector[N] c = alpha * L * z;
}

model {
  a ~ normal(0, 1.0);
  b ~ normal(0, 1.0);
  z ~ normal(0, 1.0);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ normal(a * vk0 + b * vk1 + c, sigma);
}

generated quantities {
  vector[N] yhat;
  vector[N] f = a * vk0 + b * vk1;
  
  {
    vector[N] mu = f + c;
    for(n in 1:N)
      yhat[n] = normal_rng(mu[n], sigma);
  }
}