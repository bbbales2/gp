functions {
  real[] sho(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    dydt[1] = y[2];
    dydt[2] = theta[1] * sin(y[1]);
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
  real c;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
}

transformed parameters {
  real theta[1];
  real yh[N, 2];
  
  theta[1] = c;
  
  yh = integrate_ode_rk45(sho, y0, 0.0, t, theta, x_r, x_i);
}

model {
  c ~ normal(0.0, 10.0);
  
  y ~ normal(yh[:, 1], sigmay);
  yp ~ normal(yh[:, 2], sigmayp);
}

generated quantities {
  vector[N] ymu = to_vector(yh[:, 1]);
  vector[N] ypmu = to_vector(yh[:, 2]);
}
