library(MASS)
library(tidyverse)
library(rstan)
library(ggplot2)
library(deSolve)

setwd("~/gp")

#l = 1.0
#a = 1.0
#s = 1e-1

gl = 9.8

source("derivative_kernels.R")

full = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * sin(state[['y']])))
}

linearized = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * state[['y']]))
}

times = seq(0, 10, by = 0.05)
y0 = c(y = pi / 4.0, yp = 0.0)
fout = as_tibble(ode(y = y0, times = times, func = full, parms = c(gl = gl))[,])
lout = as_tibble(ode(y = y0, times = times, func = linearized, parms = c(gl = gl))[,])

df = bind_rows(lout %>% mutate(type = "linearized"), fout %>% mutate(type = "full")) %>% mutate(type = factor(type))

df = df %>% mutate(ynoise = y + rnorm(nrow(df), mean = 0.0, sd = 0.1)) %>%
  mutate(ypnoise = yp + rnorm(nrow(df), mean = 0.0, sd = 0.1))

#df = df %>% mutate(ypl = -gl * y)

dft = df %>% filter(type == "full")
#%>%
# group_by(type) %>% mutate(ypfd = (lead(ynoise) - lag(ynoise)) / (lead(time) - lag(time))) %>% ungroup()

dft %>% ggplot(aes(time, yp)) +
  geom_line() +
  geom_point(aes(time, ypfd)) +
  geom_point(aes(time, ypl))

df %>% ggplot(aes(time, y)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(time, ynoise, colour = type))

sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise)

fit = stan("/home/bbales2/gp/models/fit_hyperparameters.stan", data = sdata, cores = 4, iter = 1000)

s = extract(fit, pars = c("rho", "alpha", "sigma"))

pairs(s)

l = median(s$rho)
a = median(s$alpha)
sy = median(s$sigma)

N = nrow(dft)
ti = dft$time
y = c(dft$ynoise, dft$ypnoise)
 
KK = matrix(0, N, N)
K = matrix(0, 2 * N, 2 * N)
KsK = matrix(0, 2 * N, 2 * N)
KsKs = matrix(0, 2 * N, 2 * N)

for(jt in 1:N) {
  for(kt in 1:N) {
    K[jt, kt] = QQ(ti[jt], ti[kt])
    K[jt + N, kt] = RQ(ti[jt], ti[kt])
    K[jt, kt + N] = QR(ti[jt], ti[kt])
    K[jt + N, kt + N] = RR(ti[jt], ti[kt])
    
    KsK[jt, kt] = RQ(ti[jt], ti[kt])
    KsK[jt + N, kt] = TQ(ti[jt], ti[kt])
    KsK[jt, kt + N] = RR(ti[jt], ti[kt])
    KsK[jt + N, kt + N] = TR(ti[jt], ti[kt])
    
    KsKs[jt, kt] = RR(ti[jt], ti[kt])
    KsKs[jt + N, kt] = TR(ti[jt], ti[kt])
    KsKs[jt, kt + N] = RT(ti[jt], ti[kt])
    KsKs[jt + N, kt + N] = TT(ti[jt], ti[kt])
  }
}

build_mu = function(K, KsK, y) {
  diag(K) = diag(K) + sy^2
  
  KsK %*% solve(K, y)
}

build_cov = function(K, KsK, KsKs) {
  KKs = t(KsK)
  
  diag(K) = diag(K) + sy^2

  KsKs - KsK %*% solve(K, KKs) + diag(1e-8, nrow(K), ncol(K))
}

out = mvrnorm(5, build_mu(K, KsK, y), build_cov(K, KsK, KsKs))
LL = build_cov(K, KsK, KsKs)
#out1 = array(out, dim = c(5, 201, 2))
#plot(ti, out1[1,, 1])
#points(ti, out1[2,, 1])
#points(ti, out1[3,, 1])
#points(ti, out1[4,, 1])
#points(ti, y[1:201], col = "red")

sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise,
             yp = dft$ypnoise,
             mu = as.vector(build_mu(K, KsK, c(dft$ynoise, dft$ypnoise))),
             L = chol(build_cov(K, KsK, KsKs)))

fit = stan("/home/bbales2/gp/models/fit_gp.stan", data = sdata, chains = 1, iter = 200)

s1 = extract(fit, c("a", "b", "c", "d", "sigmay", "sigmayp"))

pairs(s1)

sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise,
             yp = dft$ypnoise)

fit = stan("/home/bbales2/gp/models/fit_fd.stan", data = sdata, cores = 4, iter = 200)

s1 = extract(fit, c("a", "b", "c", "d", "sigmay", "sigmayp"))

pairs(s1)

fit = stan("/home/bbales2/gp/models/fit_fd_correct.stan", data = sdata, cores = 4, iter = 1000)

s2 = extract(fit, c("b", "c", "sigmay", "sigmayp"))

pairs(s2)
