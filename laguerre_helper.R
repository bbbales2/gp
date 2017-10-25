library(orthopolynom)
library(tidyverse)
library(ggplot2)
library(rstan)

laguerre_integrate = function(m, alpha) {
  nn = m + 1
  r = glaguerre.recurrences(nn + 1, alpha = alpha, normalized = TRUE)
  m = monic.polynomial.recurrences(r)
  roots = polynomial.roots(m)
  p = glaguerre.polynomials(nn + 1, alpha = alpha, normalized = TRUE)
  w = roots[[nn]] / (nn^2 * polynomial.values(p[nn + 1], roots[[nn]])[[1]]^2)
  list(roots = roots[[nn]], weights = w, polynomial = p[[nn]])
}
