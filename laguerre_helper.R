library(orthopolynom)
library(tidyverse)
library(ggplot2)
library(rstan)

laguerre_integrate = function(m) {
  nn = m + 1
  r = laguerre.recurrences(nn + 1)
  m = monic.polynomial.recurrences(r)
  roots = polynomial.roots(m)
  p = laguerre.polynomials(nn + 1)
  w = roots[[nn]] / (nn^2 * polynomial.values(p[nn + 1], roots[[nn]])[[1]]^2)
  list(roots = roots[[nn]], weights = w, polynomial = p[[nn]])
}
