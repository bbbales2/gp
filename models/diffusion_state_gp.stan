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
  
  real[] sho(real t, real[] u, real[] siglz, real[] dxscale, int[] x_i) {
    int N = size(u);
    int M = size(siglz) - 2;
    real dx = dxscale[1];
    real scale = dxscale[2];
    real dudt[N];
    real sigmagp = siglz[1];
    real l = siglz[2];
    real uD[N + 1];
    vector[N + 1] D;
    
    uD[1] = (1.0 + u[1]) / 2.0 - 0.5;
    uD[N + 1] = (u[N] + 1.0) / 2.0 - 0.5;
    for(n in 2:N)
      uD[n] = (u[n - 1] + u[n]) / 2 - 0.5;
  
    D = exp(approx_L(M, scale, uD, sigmagp, l) * to_vector(siglz[3:]));
    
    dudt[1] = (D[2] * (u[2] - u[1]) - D[1] * (u[1] - 1.0)) / dx^2;
    dudt[N] = (D[N + 1] * (1.0 - u[N]) - D[N] * (u[N] - u[N - 1])) / dx^2;
    
    for(i in 2 : (N - 1)) {
      dudt[i] = (D[i + 1] * (u[i + 1] - u[i]) - D[i] * (u[i] - u[i - 1])) / dx^2;
    }
    
    return dudt;
  }
}

data {
  int<lower=1> N;
  int M;
  int K;
  real t;
  real dx;
  real scale;
  real u0[N];
  real uDp[K];
  vector[N] u;
}

transformed data {
  real x_r[2];
  real ts[1];
  int x_i[0];
  
  ts[1] = t;
  x_r[1] = dx;
  x_r[2] = scale;
}

parameters {
  real z[M];
  real<lower=0> sigmagp;
  real<lower=0> l;
  real<lower=0> sigma;
}

transformed parameters {
  real uh[1, N];
  
  {
    real siglz[M + 2];
    siglz[1] = sigmagp;
    siglz[2] = l;
    siglz[3:] = z;
    uh = integrate_ode_rk45(sho, u0, 0.0, ts, siglz, x_r, x_i);
  }
}

model {
  vector[6] Y[5];
  vector[6] X[5];
  matrix[6, 6] Z;
  
  Y ~ multi_normal(X, exp(Z) / 5.0);
  
  z ~ normal(0.0, 1.0);
  l ~ gamma(4.0, 4.0);
  sigmagp ~ normal(0, 1.0);
  sigma ~ normal(0, 1.0);
  
  u ~ normal(uh[1, :], sigma);
}

generated quantities {
  vector[K] D;

  {
    real uDpt[K];
    
    for(n in 1:K)
      uDpt[n] = uDp[n] - 0.5;

    D = exp(approx_L(M, scale, uDpt, sigmagp, l) * to_vector(z));
  }
}
