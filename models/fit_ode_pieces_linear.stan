
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
  real zz[N, 2];
}

transformed parameters {
  real z[N, 2];
  matrix[2, 2] A;
  
  A[1, 1] = 0.0;
  A[1, 2] = 1.0;
  A[2, 1] = c;
  A[2, 2] = 0;
  
  {
    vector[2] out = matrix_exp(t[1] * A) * to_vector(y0);
    //integrate_ode_rk45(sho, y0, 0.0, { t[1] }, theta, x_r, x_i);
    for(i in 1:2)
      z[1, i] = out[i];
    
    for(n in 2:N) {
      out = matrix_exp((t[n] - t[n - 1]) * A) * to_vector(z[n - 1]);
      //out = integrate_ode_rk45(sho, z[n - 1], t[n - 1], { t[n] }, theta, x_r, x_i);
      for(i in 1:2)
        z[n, i] = out[i];
    }
  }
}

model {
  sigmay ~ normal(0.0, 1.0);
  sigmayp ~ normal(0.0, 1.0);
  
  //c ~ normal(0.0, 10.0);
  
  y ~ normal(z[:, 1], sigmay);
  yp ~ normal(z[:, 2], sigmayp);
}

generated quantities {
  vector[N] ymu = to_vector(z[:, 1]);
  vector[N] ypmu = to_vector(z[:, 2]);
}
