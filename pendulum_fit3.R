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

df %>% ggplot(aes(time, y)) +
  geom_line(aes(colour = type)) +
  geom_point(aes(time, ynoise, colour = type), size = 0.5)

dft = df %>%
  filter(type == "full")

# Fit data using finite difference derivative approximations and linear system
sdata = list(N = nrow(dft),
             t = dft$time,
             y0 = y0,
             y = dft$ynoise,
             yp = dft$ypnoise)

fit1 = stan("models/fit_ode_pieces_linear_latent.stan", data = sdata, cores = 4, chains = 4, iter = 2000)

extract(fit1, c("zet"))$zet[,, 2] %>% t %>%
  as.tibble %>%
  mutate(time = dft$time) %>%
  gather(which, value, -time) %>%
  group_by(which) %>%
  mutate(value = cumsum(value)) %>%
  ggplot(aes(time, value)) +
  geom_point(aes(group = which), alpha = 0.1, size = 0.5) +
  geom_point(data = df, aes(time, ynoise, colour = type), size = 0.5)

fit2 = stan("models/fit_ode_pieces_linear.stan", data = sdata, cores = 4, chains = 4, iter = 2000)

extract(fit2, c("zz"))$zz[,, 1] %>% t %>%
  as.tibble %>%
  mutate(time = dft$time) %>%
  gather(which, value, -time) %>%
  ggplot(aes(time, value)) +
  geom_point(aes(group = which), alpha = 0.1) +
  geom_point(data = df, aes(time, ynoise, colour = type), size = 0.5)

s1 = as_tibble(extract(fit_fd, c("c")))

s1 %>% ggplot(aes(c)) +
  geom_histogram() +
  geom_vline(xintercept = -9.8, col = "darkred") +
  labs(title = "Fit finite differences using linear model")
