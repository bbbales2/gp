#%%

import pystan
import matplotlib.pyplot as plt
import numpy

#%%

N = 20
S = 2

tf = numpy.linspace(0, numpy.pi * 2.0, S * N)
t = tf[::S]
x = numpy.sin(tf)
dx = numpy.cos(tf)
dxn = dx[::S] + numpy.random.randn(N) * 0.15

plt.plot(tf, dx, '+r')
plt.plot(t, dxn, '+g')
plt.legend(['no noise', 'w noise'])
plt.show()

#%%

model = """
functions {
  real cov(real ti, real tj, real l2) {
    return exp(-(ti - tj)^2 / l2);
  }

  real covd(real ti, real tj, real l2) {
    return 2 * exp(-(ti - tj)^2 / l2) * (ti - tj) / l2;
  }

  real covdd(real ti, real tj, real l2) {
    return 2 * exp(-(ti - tj)^2 / l2) * (l2 - 2 * (ti - tj)^2) / l2^2;
  }
}

data {
  int<lower=1> N;
  int<lower=1> M;
  vector[N] t;
  vector[N] dx;
  vector[M] tp;
}

transformed data {
  vector[N] zeros;
  vector[N + 1] xdx;

  for(n in 1:N) zeros[n] = 0.0;

  xdx[1] = 0.0;
  for(n in 1:N) xdx[n + 1] = dx[n];
}

parameters {
  real<lower = 0.0> sf2;
  real<lower = 0.0> l2;
  real<lower = 0.0> s2;
}

model {
  matrix[N, N] Sigma;

  sf2 ~ cauchy(0, 5);
  l2 ~ cauchy(0, 40);
  s2 ~ cauchy(0, 5);

  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i, j] = sf2 * covdd(t[i], t[j], l2);
      Sigma[j, i] = Sigma[i, j];
    }
  }

  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = sf2 * covdd(t[k], t[k], l2) + s2;

  dx ~ multi_normal(zeros, Sigma);
}

generated quantities {
  vector[M] dxh;
  vector[M] xh;

  {
    matrix[N + 1, N + 1] Sigma;
    matrix[M, M] SigmaPP;
    matrix[M, M] SigmaIPP;
    vector[N + 1] Kinvy;
    matrix[M, N + 1] SigmaP;
    matrix[M, N + 1] SigmaIP;

    // off-diagonal elements
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma[i + 1, j + 1] = sf2 * covdd(t[i], t[j], l2);
        Sigma[j + 1, i + 1] = Sigma[i + 1, j + 1];
      }
    }

    // diagonal elements
    for (k in 1:N)
      Sigma[k + 1, k + 1] = sf2 * covdd(t[k], t[k], l2) + s2;

    for (j in 1:N) {
      Sigma[1, j + 1] = sf2 * covd(t[1], t[j], l2);
      Sigma[j + 1, 1] = Sigma[1, j + 1];
    }

    Sigma[1, 1] = sf2 * cov(0.0, 0.0, l2) + s2;

    // off-diagonal elements
    for (i in 1:(M-1)) {
      for (j in (i+1):M) {
        SigmaPP[i, j] = sf2 * covdd(tp[i], tp[j], l2);
        SigmaPP[j, i] = SigmaPP[i, j];
      }
    }

    // diagonal elements
    for (k in 1:M)
      SigmaPP[k, k] = sf2 * covdd(tp[k], tp[k], l2) + 1e-5;

    // off-diagonal elements
    for (i in 1:(M-1)) {
      for (j in (i+1):M) {
        SigmaIPP[i, j] = sf2 * cov(tp[i], tp[j], l2);
        SigmaIPP[j, i] = SigmaIPP[i, j];
      }
    }

    // diagonal elements
    for (k in 1:M)
      SigmaIPP[k, k] = sf2 * cov(tp[k], tp[k], l2) + 1e-5;

    for (i in 1:M) {
      for (j in 1:N) {
        SigmaP[i, j + 1] = sf2 * covdd(tp[i], t[j], l2);
      }

      SigmaP[i, 1] = sf2 * covd(0.0, tp[i], l2);
    }

    for (i in 1:M) {
      for (j in 1:N) {
        SigmaIP[i, j + 1] = sf2 * covd(tp[i], t[j], l2);
      }

      SigmaIP[i, 1] = sf2 * cov(tp[i], 0.0, l2);
    }

    Kinvy = Sigma \ xdx;

    dxh = multi_normal_rng(SigmaP * Kinvy, SigmaPP - SigmaP * (Sigma \ SigmaP'));
    xh = multi_normal_rng(SigmaIP * Kinvy, SigmaIPP - SigmaIP * (Sigma \ SigmaIP'));
  }
}
"""

sm = pystan.StanModel(model_code = model)

#%%

fit = sm.sampling(data = { "N" : N,
                           "M" : len(tf),
                           "t" : t,
                           "tp" : tf,
                           "dx" : dxn })

#%%

print fit

dxh = fit.extract()['dxh']
xh = fit.extract()['xh']

plt.plot(tf, dx, 'r')
plt.plot(t, dxn, 'og')
for i in range(10):
    plt.plot(tf, dxh[-i - 1, :], '+')

plt.show()
plt.plot(tf, dx, 'r')
plt.plot(t, dxn, 'og')
for i in range(10):
    plt.plot(tf, xh[-i - 1, :], '+')