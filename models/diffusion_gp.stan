functions {
  real[] sho(real t, real[] u, real[] D, real[] dx, int[] x_i) {
    int N = size(u);
    real dudt[N];

    dudt[1] = (D[2] * (u[2] - u[1]) - D[1] * (u[1] - 1.0)) / dx[1]^2;
    dudt[N] = (D[N + 1] * (1.0 - u[N]) - D[N] * (u[N] - u[N - 1])) / dx[1]^2;

    for(i in 2 : (N - 1)) {
      dudt[i] = (D[i + 1] * (u[i + 1] - u[i]) - D[i] * (u[i] - u[i - 1])) / dx[1]^2;
    }

    return dudt;
  }
  
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
  int<lower=1> N;
  int M;
  real t;
  real scale;
  real u0[N];
  real xD[N + 1];
  vector[N] u;
}

transformed data {
  real x_r[1];
  real ts[1];
  int x_i[0];
  
  ts[1] = t;
  x_r[1] = xD[2] - xD[1];
}

parameters {
  vector[N + 1] z;
  real<lower=0> sigmagp;
  real<lower=0> l;
  real<lower=0> sigma;
}

transformed parameters {
  real uh[1, N];
  vector[N + 1] D;// = inv_logit(approx_L(M, scale, xD, sigmagp, l) * z);
  
  {
    matrix[N + 1, N + 1] Sigma = cov_exp_quad(xD, sigmagp, l);
    matrix[N + 1, N + 1] L;
    
    for(n in 1:(N + 1))
      Sigma[n, n] = Sigma[n, n] + 1e-12;
      
    L = cholesky_decompose(Sigma);
    
    D = inv_logit(L * z);
  }

  uh = integrate_ode_rk45(sho, u0, 0.0, ts, to_array_1d(D), x_r, x_i);
}

model {
  z ~ normal(0.0, 1.0);
  l ~ gamma(4.0, 4.0);
  sigmagp ~ normal(0, 1.0);
  sigma ~ normal(0, 1.0);
  
  u ~ normal(uh[1, :], sigma);
}
