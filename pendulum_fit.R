library(MASS)
library(tidyverse)
library(rstan)
library(ggplot2)
library(deSolve)
library(GGally)

setwd("~/gp")

gl = 9.8

source("derivative_kernels.R")

full = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * sin(state[['y']])))
}

linearized = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * state[['y']]))
}

h = 0.05
times = seq(0, 10, by = h)
y0 = c(y = 1.0 * pi / 4.0, yp = 0.0)
fout = as_tibble(ode(y = y0, times = times, func = full, parms = c(gl = gl))[,])
lout = as_tibble(ode(y = y0, times = times, func = linearized, parms = c(gl = gl))[,])

df = bind_rows(lout %>% mutate(type = "linearized"), fout %>% mutate(type = "full")) %>% mutate(type = factor(type))

df = df %>% mutate(ynoise = y + rnorm(nrow(df), mean = 0.0, sd = 0.01)) %>%
  mutate(ypnoise = yp + rnorm(nrow(df), mean = 0.0, sd = 0.01))

dft = df %>% filter(type == "full") %>%
  mutate(yd = (lead(y) - lag(y)) / (2.0 * h)) %>%
  mutate(ypd = (lead(yp) - lag(yp)) / (2.0 * h)) %>%
  drop_na()

df %>% ggplot(aes(time, y)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(time, ynoise, colour = type))

# Find GP hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise)
  
  fit = stan("models/fit_hyperparameters.stan", data = sdata, cores = 4, iter = 1000)
  
  sy = extract(fit, pars = c("rho", "alpha", "sigma"))
  
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ypnoise)
  
  fit = stan("models/fit_hyperparameters.stan", data = sdata, cores = 4, iter = 1000)
  
  syp = extract(fit, pars = c("rho", "alpha", "sigma"))
}

# Impute data and fit with linearized system
{
  sample_derivs = function(l, a, sy) {
    N = nrow(dft)
    ti = dft$time
    y = c(dft$ynoise, dft$ypnoise)
    
    K = matrix(0, 2 * N, 2 * N)
    KsK = matrix(0, 2 * N, 2 * N)
    KsKs = matrix(0, 2 * N, 2 * N)
    
    for(jt in 1:N) {
      for(kt in 1:N) {
        K[jt, kt] = a^2 * QQ(ti[jt], ti[kt], l)
        K[jt + N, kt] = a^2 * RQ(ti[jt], ti[kt], l)
        K[jt, kt + N] = a^2 * QR(ti[jt], ti[kt], l)
        K[jt + N, kt + N] = a^2 * RR(ti[jt], ti[kt], l)
        
        KsK[jt, kt] = a^2 * RQ(ti[jt], ti[kt], l)
        KsK[jt + N, kt] = a^2 * TQ(ti[jt], ti[kt], l)
        KsK[jt, kt + N] = a^2 * RR(ti[jt], ti[kt], l)
        KsK[jt + N, kt + N] = a^2 * TR(ti[jt], ti[kt], l)
        
        KsKs[jt, kt] = a^2 * RR(ti[jt], ti[kt], l)
        KsKs[jt + N, kt] = a^2 * TR(ti[jt], ti[kt], l)
        KsKs[jt, kt + N] = a^2 * RT(ti[jt], ti[kt], l)
        KsKs[jt + N, kt + N] = a^2 * TT(ti[jt], ti[kt], l)
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
    
    out = mvrnorm(1, build_mu(K, KsK, y), build_cov(K, KsK, KsKs))
    list(out[1:N], out[(N + 1) : (2 * N)])
  }
  
  fits = list()
  sf = stan_model("models/fit_ypypp.stan")
  for(i in 1:500) {
    ypypp = sample_derivs(s$rho[i], s$alpha[i], s$sigma[i])
  
    sdata = list(N = nrow(dft),
                 t = dft$time,
                 y = dft$ynoise,
                 yp = dft$ypnoise,
                 yph = ypypp[[1]],
                 ypph = ypypp[[2]])
    
    fits[i] = sampling(sf, data = sdata, chains = 1, cores = 1, iter = 1000)
  }
  
  df_impute2 = bind_rows(lapply(fits, function(x) {
    extract(x, pars = c("a", "b", "c", "d", "sigmay", "sigmayp"))
  })) %>% sample_frac(1.0)
  
  df_impute2 %>% sample_n(2000) %>% ggpairs()
}

# Fit data using finite difference derivative approximations and linear system
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd)
  
  fit_fd = stan("models/fit_fd.stan", data = sdata, cores = 4, iter = 2000)
  
  df_fd = as_tibble(extract(fit_fd, c("a", "b", "c", "d")))#, "sigmay", "sigmayp"
  
  df_fd %>% ggpairs
}

# Fit data using linearized ode
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               y0 = y0)
  
  fit = stan("models/fit_ode.stan", data = sdata, chains = 1, cores = 1, iter = 1000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d", "sigmayp", "sigmaypp")))
  
  s1 %>% ggpairs
}

# Fit data using linearized ode
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               y0 = y0)
  
  fit = stan("models/fit_ode.stan", data = sdata, chains = 4, cores = 4, iter = 2000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d", "sigmayp", "sigmaypp")))
  
  s1 %>% ggpairs
}

# Fit data using finite difference derivatives and Kennedy-O'Hagan-like model
#   and infer all the hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd)
  
  fit = stan("models/fit_gp.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d", "alphayp", "alphaypp", "sigmayp", "sigmaypp")))
  
  s1 %>% ggpairs
}

# Fit data using finite difference derivatives and Kennedy-O'Hagan-like model
#   and do *not* infer all the hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd,
               alphayp = 0.25,
               alphaypp = 0.25,
               rho = 2.0)
  
  fit = stan("models/fit_gp_fixed_hyperparam.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d")))#, "sigmayp", "sigmaypp"
  
  s1 %>% ggpairs
}

s1$type = 'gp'
df_fd$type = 'lr'
bind_rows(s1, df_fd) %>% ggplot(aes(c)) +
  geom_histogram(aes(fill = type), binwidth = 0.01)

get_lines = function(fit, x, vnames, n = 100) {
  a = extract(fit, vnames)
  idxs = sample(nrow(a[[vnames[[1]]]]), n)
  
  out = as_tibble()
  for(i in 1:length(vnames)) {
    vname = vnames[[i]];
    d = a[[vname]][idxs,]
    colnames(d) <- x
    d = as_tibble(d) %>% gather(time, data)
    d$name = vname
    
    out = bind_rows(out, d)
  }
  out %>% mutate(time = as.double(time))
}
a = get_lines(fit, dft$time, c("lyp", "ypm"), 50)
b = get_lines(fit, dft$time, c("lypp", "yppm"), 50)
a$state = "yd"
b$state = "ypd"
a = bind_rows(a, b)

a %>% ggplot(aes(time, data)) +
  geom_jitter(aes(color = name), width = 0.05, alpha = 0.25, size = 0.1) +
  geom_line(data = dft %>% gather(state, data, c(yd, ypd)), aes(color = state)) +
  facet_grid(state ~ .) +
  guides(colour = guide_legend(override.aes = list(size = 1, alpha = 1.0)))

# Fit data using finite difference derivatives and full non-linear model
{
  fit = stan("/home/bbales2/gp/models/fit_fd_correct.stan", data = sdata, cores = 4, iter = 1000)
  
  s2 = extract(fit, c("b", "c", "sigmay", "sigmayp"))
  
  pairs(s2)
}
