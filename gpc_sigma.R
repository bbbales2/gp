library(Rcpp)
library(tidyverse)
library(ggplot2)
library(rstan)

source("laguerre_helper.R")
sourceCpp("covariance.cpp")

M = 11 # Order of GPC expansion
K = 15 # Order of polynomial expansions for integration
k = 2.0
alpha = k - 1.0

x = seq(0.0, 1.0, length = 3)
func = function(l) {
  sin(2 * l)#c(chol(rbf_cov(x, x, l)))
}

quad = laguerre_integrate(K, alpha)

polys = list()
lookupt = list()
for(m in 0:M) {
  polys[[m + 1]] = laguerre_integrate(m, alpha)$polynomial
  output = map(quad$roots, func)
  for(i in 1:length(output)) {
    output[[i]] = output[[i]] *
      quad$weights[[i]] *
      polynomial.values(polys[m + 1], quad$roots[[i]])[[1]]
  }

  lookupt[[m + 1]] = Reduce('+', output)
}

rbf_cov_approx = function(xl) {
  Reduce('+', lapply(0:M, function(m) { lookupt[[m + 1]] * polynomial.values(polys[m + 1], xl)[[1]] }))
}

list(lvs = rgamma(100, k)) %>% as.tibble %>%
  mutate(approx = map(c(lvs), rbf_cov_approx),
         exact = map(c(lvs), func)) %>% unnest() %>%
  gather(which, y, c(-lvs)) %>%
  ggplot(aes(lvs, y)) +
  geom_point(aes(colour = which))

list(lvs = seq(0.01, 5.0, length = 100))
errors = map(lvs, function(xl) max(abs(rbf_cov_approx(xl) - func(xl))))
plot(lvs, errors)
