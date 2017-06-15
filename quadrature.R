library(orthopolynom)
library(tidyverse)
library(ggplot2)
library(rstan)

nn = 50
r = hermite.h.recurrences(nn + 1, normalized=FALSE)
m = monic.polynomial.recurrences(r)
roots = polynomial.roots(m)
p = hermite.h.polynomials(nn + 1, normalized=FALSE)
w = 2^(nn - 1) * factorial(nn) * sqrt(pi) / (nn^2 * polynomial.values(p[nn], roots[[nn + 1]])[[1]]^2)

setwd("~/gp")

source("spectralGP.R")

integrateL = function(scale, M, x, sigma, l) {
  a = 1 / (4 * scale^2)
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
  Ht[, 1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta)# * exp(-delta^2 * x * x)
  Ht[, 2] = f * sqrt(1 / 2) * 2 * xp * Ht[, 1]
  for(n in 3:M) {
    Ht[, n] = f * sqrt(1 / (2 * (n - 1))) * 2 * xp * Ht[, n - 1] - f^2 * sqrt(1 / (4 * (n - 1) * (n - 2))) * 2 * (n - 2) * Ht[, n - 2]
  }
  
  sigma * Ht
}

relu = function(x) {
  (x > 0.0) * x
}

M = 5
scale = 1.0
sigma = 1.5
l = 0.5
z = rnorm(M)
x = seq(-10.0, 10.0, length = 1000)
sum((approxL(scale, M, x, sigma, l) %*% z)) * (x[2] - x[1])
plot(x, (approxL(scale, M, x, sigma, l) %*% z))
sum(w * (approxL(scale, M, roots[[nn + 1]], sigma, l) %*% z) * exp(roots[[nn + 1]]^2))
sum(w * (integrateL(scale, M, roots[[nn + 1]] / getd(scale, sigma, l), sigma, l) %*% z)) / getd(scale, sigma, l)

a = approxL(1.0, 10, x, 1.5, 0.5)
plot(a[,4])
#colnames(a) = 1:10
#as_tibble(a) %>% gather(key, value, 1:10) %>% ggplot(aes(
