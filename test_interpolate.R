library(tidyverse)
library(ggplot2)
library(rstan)
library(Rcpp)
library(shinystan)

sourceCpp("covariance.cpp")

N = 150
x = seq(0.0, 10.0, length = N)
y = sin(x) + rnorm(N, 0.0, 0.2)
P = 10
lp = seq(qgamma(0.05, 4.0, 4.0), qgamma(0.95, 4.0, 4.0), length = 10)

Ls = list()
dLdls = list()
for(p in 1:P) {
  out = rbf_cov_chol(x, lp[p])
  Ls[[p]] = out$L
  dLdls[[p]] = out$dLdl
}

list(x = x, y = y) %>%
  as.tibble %>%
  ggplot(aes(x, y)) +
  geom_point()

cpp_src = stanc("models/cubic_interpolated_gp.stan", allow_undefined = TRUE)$cppcode
cpp_src_split = strsplit(cpp_src, "\n")[[1]]
first_match = grep("approx_Lz", cpp_src_split)[[1]]
cat(cpp_src_split[(first_match - 2) : (first_match + 4)], sep = "\n")

model = stan_model("models/cubic_interpolated_gp.stan",
                   allow_undefined = TRUE,
                   includes = paste0('\n#include "',
                                     file.path(getwd(),
                                               'models/cubic_interpolated_gp.hpp'), '"\n'))

(fit1 = sampling(model, data = list(N = N,
                                   x = x,
                                   y = y,
                                   P = P,
                                   lp = lp,
                                   Ls = Ls,
                                   dLdls = dLdls),
                iter = 1000,
                chains = 1))

(fit2 = stan("models/exact_gp.stan",
             data = list(N = N,
                         x = x,
                         y = y),
             iter = 1000,
             chains = 1))

(fit3 = stan("models/interpolated_gp.stan",
             data = list(N = N,
                         x = x,
                         y = y,
                         P = P,
                         lp, lp),
             iter = 1000,
             chains = 1))
