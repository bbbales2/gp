library(MASS)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(deSolve)
library(GGally)
library(purrr)
library(parallel)

setwd("~/gp")

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

dft = df %>% filter(type == "full")

df %>% ggplot(aes(time, y)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(time, ynoise, colour = type), size = 0.5)

####################################################
####################################################
# For a process observed with noise at a finite number of points,
# Find the hyperparameters of the SEXP kernel describing it.
# By find we mean get draws from the posterior of these hyperparameters
####################################################
####################################################
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise)
  
  fit = stan("models/fit_hyperparameters.stan", data = sdata, cores = 4, iter = 1000)
  
  s = extract(fit, pars = c("rho", "alpha", "sigma"))
}

####################################################
####################################################
# once we get posterior samples of the hyperparameters we can sample from them
# to generate posterior predictive draws of the derivative of the original process.
####################################################
####################################################
{
  #function that takes single posterior sample of rho, alpha, and sigma and returns
  #a single posterior draw of the derivative process (w/o noise)
  sample_derivs = function(params, ynoise, ti) {
    l <- params[1]
    a <- params[2]
    sy <- params[3]
    
    N = length(ynoise)

    #use outer to create kernel evaluated at all combinations of time points.
    #the original kernel function take in two time poitns AND an l so we use
    #kern_fixed_l to make versions of the kernel function w/ a fixed l
    kern_fixed_l <- function(kern, l) {function(tj, tk) kern(tj, tk, l)} 
    K <- a^2 * outer(ti, ti, FUN = kern_fixed_l(QQ, l))    #cov between old points and old points
    KsK <- a^2 * outer(ti, ti, FUN = kern_fixed_l(RQ, l))  #old points and new points
    KsKs <- a^2 * outer(ti, ti, FUN = kern_fixed_l(RR, l)) #new points and new points
    
    build_mu = function(K, KsK, y) {
      diag(K) = diag(K) + sy^2
      KsK %*% solve(K, ynoise)
    }
    
    build_cov = function(K, KsK, KsKs) {
      KKs = t(KsK)
      diag(K) = diag(K) + sy^2
      KsKs - KsK %*% solve(K, KKs) + diag(1e-8, nrow(K), ncol(K))
    }
    
    out = mvrnorm(1, build_mu(K, KsK, ynoise), build_cov(K, KsK, KsKs))
    return(out)
  }
  
  #get all of our posterior samples as a list of vectors c(rho, alpha, sigma)
  get_single_sample = function(i) c(rho = s$rho[i], alpha = s$alpha[i], sigma = s$sigma[i])
  s_list = map(1:500, get_single_sample)
  
  #for each posterior sample of c(rho, alpha, sigma) get a post. pred. sample of deriv. for both states
  sample_derivs_both_states <- function(params, ynoise, ypnoise, ti) {
    return(list(yp_hat = sample_derivs(params, ynoise, ti),
                ypp_hat = sample_derivs(params, ypnoise, ti)))
  }
  post_draws <- mclapply(s_list, sample_derivs_both_states, ynoise = dft$ynoise, ypnoise = dft$ypnoise, ti = dft$time, mc.cores = 2)
  
  #create a function that does a Stan fit given a list containing yp_hat and ypp_hat
  sf = stan_model("models/fit_ypypp.stan")
  
  get_stan_fit_from_impute <- function(imputed_draw) {
    sdata = list(N = nrow(dft),
                 t = dft$time,
                 y = dft$ynoise,
                 yp = dft$ypnoise,
                 yph = imputed_draw$yp_hat,
                 ypph = imputed_draw$ypp_hat)
    
    fit <- sampling(sf, data = sdata, chains = 1, cores = 1, iter = 1000)
    return(fit)
  }
  
  #apply function to get a single Stan fit for each posterior predictive draw all in one big list
  stan_fits <- mclapply(post_draws, get_stan_fit_from_impute, mc.cores = 4)
  
  #coalesce samples from these stan fits in to a single dataframe
  df_impute2 = bind_rows(lapply(stan_fits, function(x) {
    extract(x, pars = c("a", "b", "c", "d", "sigmay", "sigmayp"))
  })) %>% sample_frac(1.0)
  
  #pairs plot
  df_impute2 %>% sample_n(2000) %>% ggpairs()
}

####################################################
####################################################
# Impute data and fit with linearized system
####################################################
####################################################
{
  sample_derivs = function(l, a, sy) {
    N = nrow(dft)
    ti = dft$time
    y = c(dft$ynoise, dft$ypnoise)
    
    K = matrix(0, 2 * N, 2 * N)      #cov between old points and old points
    KsK = matrix(0, 2 * N, 2 * N)    #old points and new points
    KsKs = matrix(0, 2 * N, 2 * N)   #new points and new points
    
    #
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
  sf = stan_model("/home/bbales2/gp/models/fit_ypypp.stan")
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
               yp = dft$ypnoise)
  
  fit = stan("/home/bbales2/gp/models/fit_fd.stan", data = sdata, cores = 4, iter = 2000)
  
  df_fd = as_tibble(extract(fit, c("a", "b", "c", "d", "sigmay", "sigmayp")))
  
  df_fd %>% ggpairs
}

# Fit data using finite difference derivatives and Kennedy-O'Hagan-like model
#   and infer all the hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise)
  
  fit = stan("/home/bbales2/gp/models/fit_gp.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d", "alphayp", "alphaypp", "sigmayp", "sigmaypp")))
  
  plot(t(extract(fit, "ypm"))[[1]][400,])
  
  s1 %>% ggpairs
}

# Fit data using finite difference derivatives and Kennedy-O'Hagan-like model
#   and do *not* infer all the hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               alphayp = 0.1,
               alphaypp = 0.1,
               rho = 0.5)
  
  fit = stan("/home/bbales2/gp/models/fit_gp_fixed_hyperparam.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit, c("a", "b", "c", "d", "sigmayp", "sigmaypp")))
  
  plot(t(extract(fit, "ypm"))[[1]][1000,])
  
  s1 %>% ggpairs
}

# Fit data using finite difference derivatives and full non-linear model
{
  fit = stan("/home/bbales2/gp/models/fit_fd_correct.stan", data = sdata, cores = 4, iter = 1000)
  
  s2 = extract(fit, c("b", "c", "sigmay", "sigmayp"))
  
  pairs(s2)
}
