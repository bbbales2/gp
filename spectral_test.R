library(tidyverse)
library(rstan)
library(ggplot2)
library(tictoc)

bH = function(x, M, l) {
  sigma = 2.0
  a = 1 / (4 * sigma^2)
  b = 1 / (2 * l^2)
  c = sqrt(a^2 + 2 * a * b)
  
  epsilon = sqrt(b)
  alpha = sqrt(2 * a)
  beta = (1 + (2 * epsilon / alpha)^2)^0.25
  delta = sqrt(alpha^2 * (beta^2 - 1) / 2)
  
  N = length(x)
  Ht = matrix(0, nrow = N, ncol = M)
  xp = alpha * beta * x
  f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2))
  Ht[, 1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x * x)
  Ht[, 2] = f * sqrt(1 / 2) * 2 * xp * Ht[, 1]
  for(n in 3:M) {
    Ht[, n] = f * sqrt(1 / (2 * (n - 1))) * 2 * xp * Ht[, n - 1] - f^2 * sqrt(1 / (4 * (n - 1) * (n - 2))) * 2 * (n - 2) * Ht[, n - 2]
  }
  
  Ht
}

bL = function(x, l) {
  sigma = 1.0
  K = outer(x, x, function(xi, xj) { sigma^2 * exp(-(xi - xj)^2 / (2 * l^2)) })
  t(chol(K + diag(1e-12, length(x))))
}

N = 200
minx = -5.0
maxx = 5.0
l = 2.0

I = 10

logspace = function(x, y, N) {
  exp(seq(log(x), log(y), length = N))
}

#Rs = 1:5

gete = function(l, M, x_mag) {
  x = runif(N, -x_mag, x_mag)
  
  Lt = bH(x, M, l)
  L = bL(x, l)
  
  K = L %*% t(L)
  Kt = Lt %*% t(Lt)
  
  max(abs(K - Kt))
}

#Ls = logspace(0.25, 3.0, I)
Xs = seq(1.0, 10.0, length = 50)
Ls = seq(0.05, 2.0, length = 50)
Ms = seq(10, 90, 10)
x = runif(N, minx, maxx)

df = as_tibble(expand.grid(Ls, Ms, Xs)) %>%
  rename(l = Var1, M = Var2, x_mag = Var3)

ef = df %>% rowwise %>% mutate(error = gete(l, M, x_mag)) %>% ungroup()

ef %>% ggplot(aes(l, error)) +
  geom_line(aes(group = x_mag, color = x_mag)) +
  scale_y_log10() +
  facet_wrap(~ M)

ef %>% filter(l == Ls[[12]]) %>% ggplot(aes(M, error)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~ x_mag, labeller = "label_both")

ef %>% filter(l == Ls[[12]]) %>% ggplot(aes(x_mag, error)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~ M, labeller = "label_both")

ef %>% filter(error < 1e-10) %>%
  group_by(l, M) %>%
  summarize(max_x = max(x_mag)) %>% ungroup() %>%
  ggplot(aes(max_x, M)) +
  scale_y_log10() +
  geom_line(aes(group = l, color = l)) +
  scale_colour_gradient2()

ef %>% mutate(xl = x_mag / l^2) %>% ggplot(aes(xl, error)) +
  geom_point(aes(color = l)) +
  facet_wrap(~ M) +
  scale_x_log10() +
  scale_y_log10()

system.time(replicate(100, bH(x, 50, l)))
system.time(replicate(100, bL(x, l)))
