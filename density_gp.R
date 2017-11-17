library(tidyverse)
library(rstan)
library(ggplot2)
library(shinystan)

setwd("~/gp")

N = 20
M = 20
I = 500
scale = 0.35
x = c(rnorm(N / 2) * 0.25 - 0.5, rnorm(N / 2) * 0.25 + 0.5)
xi = seq(min(x) - 0.5, max(x) + 0.5, length = I)

fit = stan("models/density_gp.stan", data = list(N = N,
                                           M = M,
                                           I = I,
                                           scale = scale,
                                           x = x,
                                           xi = xi), cores = 4, iter = 1000)

sm = stan_model("models/density_gp.stan")

fit2 = optimizing(sm, data = list(N = N,
                                  M = M,
                                  I = I,
                                  scale = scale,
                                  x = x,
                                  xi = xi))

fit = vb(sm, data = list(N = N,
                          M = M,
                          I = I,
                          scale = scale,
                          x = x,
                          xi = xi))

a = as_tibble(extract(fit, "f")$f)
colnames(a) = 1:ncol(a)
b = as_tibble(list(idx = 1:length(x), x = x))
out = inner_join(a %>% gather(idx, f) %>% mutate(idx = as.integer(idx)), b, by = "idx")

summary = out %>% group_by(x) %>%
  summarize(mean = mean(f),
            q1 = quantile(f, 0.025),
            q2 = quantile(f, 0.167),
            q3 = quantile(f, 0.833),
            q4 = quantile(f, 0.975)) %>%
  ungroup()

summary %>% ggplot(aes(x, mean)) +
  geom_ribbon(aes(ymin = q1, ymax = q2), alpha = 0.25) +
  geom_ribbon(aes(ymin = q2, ymax = q3), alpha = 0.75) +
  geom_ribbon(aes(ymin = q3, ymax = q4), alpha = 0.25) +
  geom_line() +
  geom_line(aes(x, q1), alpha = 1.0, size = 0.125) +
  geom_line(aes(x, q2), alpha = 1.0, size = 0.125) +
  geom_line(aes(x, q3), alpha = 1.0, size = 0.125) +
  geom_line(aes(x, q4), alpha = 1.0, size = 0.125) +
  #geom_point(data = out2, aes(time, f), size = 0.1, alpha = 0.01) +
  #geom_boxplot(data = summary2, aes(time, f, group = time)) +
  xlab("x") +
  ylab("pdf") +
  ggtitle("estimated pdf (95% confidence intervals blue)") +
  geom_rug(data = as_tibble(list(x = x)), aes(x = x, y = NULL), colour = "purple") +
  stat_function(fun = function(x) {
    dnorm(x, sd = 0.25, mean = 0.5) / 2.0 + dnorm(x, sd = 0.25, mean = -0.5) / 2.0
    }, colour = "red") +
  geom_line(data = as_tibble(list(x = x, y = fit2$par[grep("f", names(fit2$par))][1:20])), aes(x, y), colour = "green")

