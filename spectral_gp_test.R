library(MASS)
library(tidyverse)
library(rstan)
library(ggplot2)
library(deSolve)
library(GGally)

setwd("~/gp")

gl = 9.8

N = 1000
x = seq(-2.0, 2.0, length = N)

alpha = 1.0;
l = 0.15;
sigma = 0.32

y = mvrnorm(1, rep(0, N), outer(x, x, function(xi, xj) { alpha^2 * exp(-(xi - xj)^2 / (2 * l^2)) + diag(sigma, N) }))

df = as_tibble(list(x = x, y = y))

Nsub = 200
idxs = sample(1:N, Nsub)
dfs = df[idxs, ]

df %>% ggplot(aes(x, y)) +
  geom_point(alpha = 0.5) +
  geom_point(data = dfs, col = "red")

# Find GP hyperparameters
sdata = list(N = nrow(dfs),
             M = 10,
             x = dfs$x,
             y = dfs$y)
  
fit = stan("models/fit_approx_gp.stan", data = sdata, chains = 1, cores = 1, iter = 1000)

fitf = stan("models/fit_full_gp.stan", data = sdata, chains = 1, cores = 1, iter = 1000)

sp = as.tibble(extract(fit, pars = c("l", "alpha")))
sf = as.tibble(extract(fitf, pars = c("l", "alpha")))

sp$which = "approx"
sf$which = "exact"

dfz = bind_rows(sp, sf) %>% gather(var, data, c(l, alpha))

dfz %>% ggplot(aes(data)) +
  geom_freqpoly(aes(color = which)) +
  facet_wrap(~ var, scales = "free")

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

a = get_lines(fit, dfs$x, c("zh"), 50)
b = get_lines(fitf, dfs$x, c("z"), 50)
a$which = "approx"
b$which = "exact"
dfp = bind_rows(a, b)

dfp %>% ggplot(aes(time, data)) +
  geom_point(aes(color = name), size = 1.0, alpha = 0.2) +
  geom_line(data = dfs, aes(x, y))
