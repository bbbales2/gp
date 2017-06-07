functions {
  matrix approx_L(int N, int M, real approx_sigma, real[] xt, real sigma, real l) {
    real a = 1.0 / (4.0 * approx_sigma^2);
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
  int<lower=1> N;
  real t[N];
  vector[N] y;
  vector[N] yp;
  vector[N] yd;
  vector[N] ypd;
}

transformed data {
  real ya[N];
  real ypa[N];
  int M = 10;
  
  for(n in 1 : N) {
    ya[n] = y[n];
    ypa[n] = yp[n];
  }
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmayp;
  real<lower=0> sigmaypp;
  
  vector<lower=0.0>[2] rho;
  real<lower=0.0> alphayp;
  real<lower=0.0> alphaypp;
  
  vector[M] etayp;
  vector[M] etaypp;
}

transformed parameters {
  vector[N] ypm;
  vector[N] yppm;
  
  {
    ypm = (approx_L(N, M, 1.0, ya, alphayp, rho[1]) .* approx_L(N, M, 1.0, ypa, alphayp, rho[2])) * etayp;
    yppm = (approx_L(N, M, 1.0, ya, alphaypp, rho[1]) .* approx_L(N, M, 1.0, ypa, alphaypp, rho[2])) * etaypp;
  }
}

model {
  rho ~ gamma(4.0, 4.0);
  alphayp ~ normal(0, 1.0);
  alphaypp ~ normal(0, 1.0);
  sigmayp ~ normal(0, 1.0);
  sigmaypp ~ normal(0, 1.0);
  
  etayp ~ normal(0.0, 1.0);
  etaypp ~ normal(0.0, 1.0);
  
  yd ~ normal(a * y + b * yp + ypm, sigmayp);
  ypd ~ normal(c * y + d * yp + yppm, sigmaypp);
}

