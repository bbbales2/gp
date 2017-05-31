functions {
  real[] sho(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    dydt[1] = theta[1] * y[1] + theta[2] * y[2];
    dydt[2] = theta[3] * y[1] + theta[4] * y[2];
    return dydt;
  }
}

data {
  int<lower=1> N;
  real t[N];
  real y0[2];
  vector[N] y;
  vector[N] yp;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real a;
  real b;
  real c;
  real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

transformed parameters {
  real theta[4];
  real yh[N, 2];
  
  theta[1] = a;
  theta[2] = b;
  theta[3] = c;
  theta[4] = d;
  
  yh = integrate_ode_rk45(sho, y0, 0.0, t, theta, x_r, x_i);
}

model {
  a ~ normal(0.0, 10.0);
  b ~ normal(0.0, 10.0);
  c ~ normal(0.0, 10.0);
  d ~ normal(0.0, 10.0);
  
  y ~ normal(yh[:, 1], sigmay);
  yp ~ normal(yh[:, 2], sigmayp);
}

generated quantities {
  vector[N] ymu = to_vector(yh[:, 1]);
  vector[N] ypmu = to_vector(yh[:, 2]);
}
