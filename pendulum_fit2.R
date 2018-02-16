library(MASS)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(deSolve)
library(GGally)
library(purrr)
library(parallel)
library(assist)
library(shinystan)

gl = 9.8

full = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * sin(state[['y']])))
}

linearized = function(t, state, parameters) {
  list(c(state[['yp']], -parameters[['gl']] * state[['y']]))
}

h = 0.05
times = seq(0, 4, by = h)
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

# Fit data using finite difference derivatives and full non-linear model
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

# Fit data using full ode
sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise,
             yp = dft$ypnoise,
             y0 = y0)
  
fit_fode = stan("models/fit_ode_full.stan", data = sdata, chains = 1, cores = 1, iter = 2000)

#launch_shinystan(fit_fode)

s1 = as_tibble(extract(fit_fode, c("c")))
  
s1 %>% ggplot(aes(c)) +
  geom_histogram() +
  geom_vline(xintercept = -gl, col = "darkred") +
  labs(title = "Fit raw data using full ODE")

# Fit data using full ode piece-wise with latent states
sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise,
             yp = dft$ypnoise,
             y0 = y0)

fit_flatent = stan("models/fit_ode_pieces_latent.stan", data = sdata,
                   chains = 4, cores = 4, iter = 1000,
                   control = list(max_treedepth = 8))

launch_shinystan(fit_flatent)

s1 = as_tibble(extract(fit_flatent, c("c")))

s1 %>% ggplot(aes(c)) +
  geom_histogram() +
  geom_vline(xintercept = -gl, col = "darkred") +
  labs(title = "Fit raw data using full ODE, piecewise, with latent states")

# Fit data using linear ode piece-wise with latent states
sdata = list(N = nrow(dft),
             t = dft$time,
             y = dft$ynoise,
             yp = dft$ypnoise,
             y0 = y0)

fit_latent = stan("models/fit_ode_pieces_linear_latent.stan", data = sdata,
                  chains = 4, cores = 4, iter = 1000,
                  control = list(max_treedepth = 8))

launch_shinystan(fit_latent)

s1 = as_tibble(extract(fit_latent, c("c")))

s1 %>% ggplot(aes(c)) +
  geom_histogram() +
  geom_vline(xintercept = -gl, col = "darkred") +
  labs(title = "Fit raw data using linear ODE, piecewise, with latent states")
