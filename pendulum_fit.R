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

h = 0.05
times = seq(0, 10, by = h)
y0 = c(y = 1.0 * pi / 4.0, yp = 0.0)
fout = as_tibble(ode(y = y0, times = times, func = full, parms = c(gl = gl))[,])
lout = as_tibble(ode(y = y0, times = times, func = linearized, parms = c(gl = gl))[,])

df = bind_rows(lout %>% mutate(type = "linearized"), fout %>% mutate(type = "full")) %>% mutate(type = factor(type))

df = df %>% mutate(ynoise = y + rnorm(nrow(df), mean = 0.0, sd = 0.1)) %>%
  mutate(ypnoise = yp + rnorm(nrow(df), mean = 0.0, sd = 0.1))

dft = df %>% filter(type == "full") %>%
  mutate(yd = (lead(y) - lag(y)) / (2.0 * h)) %>%
  mutate(ypd = (lead(yp) - lag(yp)) / (2.0 * h)) %>%
  drop_na()

df %>% ggplot(aes(time, y)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(time, ynoise, colour = type), size = 0.5)


# Fit data using finite difference derivative approximations and linear system
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd)
  
  fit_fd = stan("models/fit_fd.stan", data = sdata, cores = 4, iter = 2000)
  
  s1 = as_tibble(extract(fit_fd, c("a", "b", "c", "d")))
  
  s1 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit finite differences using linear model")
}

# Fit data using finite difference derivatives and full non-linear model
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd)
  
  fit_fdfull = stan("/home/bbales2/gp/models/fit_fd_full.stan", data = sdata, cores = 4, iter = 2000)
  
  s1 = as_tibble(extract(fit_fdfull, c("c", "sigmayp")))
  
  s1 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit finite differences using full model")
}

# Fit data using linearized ode
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               y0 = y0)
  
  fit_lode = stan("models/fit_ode.stan", data = sdata, chains = 1, cores = 1, iter = 1000)
  
  s1 = as_tibble(extract(fit_lode, c("a", "b", "c", "d", "sigmay", "sigmayp")))
  
  s1 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit raw data using linearized ODE")
}

# Fit data using full ode
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               y0 = y0)
  
  fit_fode = stan("models/fit_ode_full.stan", data = sdata, chains = 1, cores = 1, iter = 1000)
  
  s1 = as_tibble(extract(fit_fode, c("c")))
  
  s1 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit raw data using full ODE")
  
  #get_lines(fit_fode, dft$time, c("ymu"), 100) %>% ggplot(aes(time, data)) +
  #  geom_point(alpha = 0.1) +
  #  geom_line(data = dft, aes(time, y), col = "red")
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
  
  fit_ko_hyp = stan("models/fit_gp.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit_ko_hyp, c("a", "b", "c", "d", "alphayp", "alphaypp", "sigmayp", "sigmaypp")))
  
  s1 %>% filter(c < -5.0) %>% ggplot(aes(c)) +
    geom_histogram(fill = "red") +
    geom_histogram(data = s2, alpha = 0.5, fill = "green") +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit finite differences using linear model + Kennedy O'Hagan + estimate hyperparameters")
}

# Fit data using finite difference derivatives and Kennedy-O'Hagan-like model
#   with approximate GPs and infer all the hyperparameters
{
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ynoise,
               yp = dft$ypnoise,
               yd = dft$yd,
               ypd = dft$ypd)
  
  fit_ko_approx = stan("models/fit_ko_approx.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s2 = as_tibble(extract(fit_ko_approx, c("a", "b", "c", "d", "alphayp", "alphaypp", "sigmayp", "sigmaypp")))
  
  s2 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit finite differences using linear model + Kennedy O'Hagan + approx GP + estimate hyperparameters")
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
  
  fit_gp_fixed = stan("models/fit_gp_fixed_hyperparam.stan", data = sdata, chains = 4, cores = 4, iter = 1000)
  
  s1 = as_tibble(extract(fit_gp_fixed, c("a", "b", "c", "d")))
  
  s1 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit finite differences using linear model + Kennedy O'Hagan + fixed hyperparameters")
}

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
  
  sy = extract(fit, pars = c("rho", "alpha", "sigma"))
  
  sdata = list(N = nrow(dft),
               t = dft$time,
               y = dft$ypnoise)
  
  fit = stan("models/fit_hyperparameters.stan", data = sdata, cores = 4, iter = 1000)
  
  syp = extract(fit, pars = c("rho", "alpha", "sigma"))
}

####################################################
####################################################
# once we get posterior samples of the hyperparameters we can sample from them
# to generate posterior predictive draws of the derivative of the original process.
####################################################
####################################################
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

{  
  #get all of our posterior samples as a list of vectors c(rho, alpha, sigma)
  get_single_sample = function(i) c(rho_y = sy$rho[i], alpha_y = sy$alpha[i], sigma_y = sy$sigma[i],
                                    rho_yp = syp$rho[i], alpha_yp = syp$alpha[i], sigma_yp = syp$sigma[i])
  s_list = map(1:100, get_single_sample)

  #for each posterior sample of c(rho, alpha, sigma) get a post. pred. sample of deriv. for both states
  sample_derivs_both_states <- function(params, ynoise, ypnoise, ti) {
    return(list(yp_hat = sample_derivs(params[1:3], ynoise, ti),
                ypp_hat = sample_derivs(params[4:6], ypnoise, ti)))
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

# Do imputation with approximate Kennedy-O'Hagan
{  
  #get all of our posterior samples as a list of vectors c(rho, alpha, sigma)
  get_single_sample = function(i) c(rho = sy$rho[i], alpha = sy$alpha[i], sigma = sy$sigma[i])
  s_list = map(1:100, get_single_sample)
  
  #for each posterior sample of c(rho, alpha, sigma) get a post. pred. sample of deriv. for both states
  sample_derivs_both_states <- function(params, ynoise, ypnoise, ti) {
    return(list(yp_hat = sample_derivs(params, ynoise, ti),
                ypp_hat = sample_derivs(params, ypnoise, ti)))
  }
  post_draws <- mclapply(s_list, sample_derivs_both_states, ynoise = dft$ynoise, ypnoise = dft$ypnoise, ti = dft$time, mc.cores = 2)
  
  #create a function that does a Stan fit given a list containing yp_hat and ypp_hat
  sf = stan_model("models/fit_ko_approx.stan")
  
  get_stan_fit_from_impute <- function(imputed_draw) {
    sdata = list(N = nrow(dft),
                 t = dft$time,
                 y = dft$ynoise,
                 yp = dft$ypnoise,
                 yd = imputed_draw$yp_hat,
                 ypd = imputed_draw$ypp_hat)
    
    fit <- sampling(sf, data = sdata, chains = 1, cores = 1, iter = 1000)
    return(fit)
  }
  
  #apply function to get a single Stan fit for each posterior predictive draw all in one big list
  stan_fits2 <- mclapply(post_draws, get_stan_fit_from_impute, mc.cores = 4)
  
  #coalesce samples from these stan fits in to a single dataframe
  df_impute3 = bind_rows(lapply(stan_fits2, function(x) {
    extract(x, pars = c("a", "b", "c", "d", "sigmayp", "sigmaypp"))
  })) %>% sample_frac(1.0)
  
  #pairs plot
  df_impute3 %>% sample_n(2000) %>% ggpairs()
  
  df_impute3 %>% ggplot(aes(c)) +
    geom_histogram() +
    geom_vline(xintercept = -9.8, col = "darkred") +
    labs(title = "Fit inputed derivatives using linear model + Kennedy O'Hagan + approx GP + estimate hyperparameters")
}

bind_rows(df_impute3 %>% sample_n(2000) %>% select(c) %>% mutate(name = "impute_ko_approx"),
          df_impute2 %>% sample_n(2000) %>% select(c) %>% mutate(name = "impute_linear"),
          as_tibble(extract(fit_ko_approx, c("c"))) %>% mutate(name = "fd_ko_approx"),
          as_tibble(extract(fit_gp_fixed, c("c"))) %>% mutate(name = "fd_ko_fixed_hp"),
          as_tibble(extract(fit_lode, c("c"))) %>% mutate(name = "ode_linear"),
          as_tibble(extract(fit_fdfull, c("c"))) %>% mutate(name = "fd_full"),
          as_tibble(extract(fit_fd, c("c"))) %>% mutate(name = "fd_linear")) %>% ggplot(aes(c)) +
  geom_histogram(aes(y = ..density..)) +
  geom_vline(xintercept = -9.8) +
  facet_grid(name ~ .)

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
