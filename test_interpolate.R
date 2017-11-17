library(tidyverse)
library(ggplot2)
library(rstan)

N = 100
x = seq(0.0, 10.0, length = N)
y = sin(x) + rnorm(N, 0.0, 0.2)
P = 10
lp = seq(qgamma(0.05, 4.0, 4.0), qgamma(0.95, 4.0, 4.0), length = 10)

list(x = x, y = y) %>%
  as.tibble %>%
  ggplot(aes(x, y)) +
  geom_point()

cpp_src = stanc("models/linear_interpolated_gp.stan", allow_undefined = TRUE)$cppcode
cpp_src_split = strsplit(cpp_src, "\n")[[1]]
first_match = grep("approx_L", cpp_src_split)[[1]]
cat(cpp_src_split[(first_match - 2) : (first_match + 3)], sep = "\n")

model = stan_model("models/linear_interpolated_gp.stan",
           allow_undefined = TRUE,
           includes = paste0('\n#include "', file.path(getwd(), 'models/linear_interpolated_gp.hpp'), '"\n'))

optimizing(model, data = list(N = N,
                              x = x,
                              y = y,
                              P = P,
                              lp, lp))

(fit1 = stan("models/exact_gp.stan",
           data = list(N = N,
                       x = x,
                       y = y),
           iter = 1000,
           chains = 1))

(fit2 = stan("models/interpolated_gp.stan",
            data = list(N = N,
                        x = x,
                        y = y,
                        P = P,
                        lp, lp),
            iter = 1000,
            chains = 1))
