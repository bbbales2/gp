library(tidyverse)
library(ggplot2)
library(rstan)

N = 150
x = seq(0.0, 10.0, length = N)
y = sin(x) + rnorm(N, 0.0, 0.2)
P = 20
lp = rgamma(P, 4.0, 4.0)

list(x = x, y = y) %>%
  as.tibble %>%
  ggplot(aes(x, y)) +
  geom_point()

(fit1 = stan("models/exact_gp.stan",
           data = list(N = N,
                       x = x,
                       y = y),
           iter = 400,
           chains = 1))

(fit2 = stan("models/interpolated_gp.stan",
            data = list(N = N,
                        x = x,
                        y = y,
                        P = P,
                        lp, lp),
            iter = 400,
            chains = 1))
