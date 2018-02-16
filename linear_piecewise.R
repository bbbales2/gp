library(tidyverse)
library(ggplot2)
library(shinystan)
library(rstan)

gl = 1.7

t = seq(0.1, 1.0, length = 10)
y = t * gl

t = c(t, t + max(t))
y = c(y, y + max(y) + 2.0)

y = y + rnorm(length(y), 0.0, 0.1)

plot(t, y)

# Fit data using finite difference derivative approximations and linear system
sdata = list(N = length(t),
             t = t,
             y = y)

fit = stan("models/linear_piecewise.stan", data = sdata, chains = 1, cores = 1, iter = 1000)

s1 = as_tibble(extract(fit, c("c")))

s1 %>% ggplot(aes(c)) +
  geom_histogram() +
  geom_vline(xintercept = gl, col = "darkred") +
  labs(title = "piecewise linear")

(extract(fit, c("z"))$z) %>%
  as.tibble %>% sample_frac(0.1) %>% setNames(t) %>%
  gather(time, value) %>%
  mutate(time = as.numeric(time)) %>%
  ggplot(aes(time, value)) +
  geom_point(alpha = 0.1) +
  geom_point(data = sdata %>% as.tibble, aes(t, y), color = "red")

launch_shinystan(fit)  
