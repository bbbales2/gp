data {
  int<lower=1> N;
  real t[N];
  vector[N] y;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real c;
  real<lower=0> sigma;
  real tau;
  vector<lower=0>[N] epsilonz;
  vector[N] zz;
}

transformed parameters {
  vector[N] zet = tau * zz .* epsilonz;
  vector[N] z;
  //vector[N] epsilon = epsilonz;
  
  z[1] = c * t[1] + zet[1];

  for(n in 2:N) {
    z[n] = z[n - 1] + c * (t[n] - t[n - 1]) + zet[n];
  }
}

model {
  tau ~ normal(0.0, 0.1);
  epsilonz ~ cauchy(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  
  c ~ normal(0.0, 10.0);
  
  zz ~ normal(0, 1.0);
  y ~ normal(z, sigma);
}
