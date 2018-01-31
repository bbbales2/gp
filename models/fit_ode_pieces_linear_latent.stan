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
  //real a;
  //real b;
  real c;
  //real d;
  real<lower=0> sigmay;
  real<lower=0> sigmayp;
  real<lower=0> epsilon;
  real z[N, 2];
}

transformed parameters {
  real theta[4];
  real zmu[N, 2];
  
  
  theta[1] = 0;
  theta[2] = 1.0;
  theta[3] = c;
  theta[4] = 0;
  
  {
    real out[1, 2] = integrate_ode_rk45(sho, y0, 0.0, { t[1] }, theta, x_r, x_i);
    for(i in 1:2)
      zmu[1, i] = out[1, i];
    
    for(n in 2:N) {
      out = integrate_ode_rk45(sho, z[n - 1], t[n - 1], { t[n] }, theta, x_r, x_i);
      for(i in 1:2)
        zmu[n, i] = out[1, i];
    }
  }
}

model {
  to_array_1d(z) ~ normal(to_array_1d(zmu), epsilon);
  
  epsilon ~ normal(0.0, 1.0);
  sigmay ~ normal(0.0, 1.0);
  sigmayp ~ normal(0.0, 1.0);
  
  //a ~ normal(0.0, 10.0);
  //b ~ normal(0.0, 10.0);
  c ~ normal(0.0, 10.0);
  //d ~ normal(0.0, 10.0);
  
  y ~ normal(z[:, 1], sigmay);
  yp ~ normal(z[:, 2], sigmayp);
}

generated quantities {
  vector[N] ymu = to_vector(z[:, 1]);
  vector[N] ypmu = to_vector(z[:, 2]);
}
