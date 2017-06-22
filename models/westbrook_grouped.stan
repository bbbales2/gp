data {
  int N;
  int T;
  int x[N];
  int y[N];
}

parameters {
  vector[T] q;
}

model {
  q ~ normal(0, 1);
  
  for(n in 1:N)
    y[n] ~ bernoulli_logit(q[x[n]]);
}

generated quantities {
  vector[T] f = inv_logit(q);
}