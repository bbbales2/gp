library(MASS)
library(reshape2)
library(tidyverse)
library(ggplot)
library(rstan)

l = 0.5
sigmaf = 1.0
sigman = 0.01
N = 10


Sigma = matrix(0, nrow = N, ncol = N)

for(i in 1:N) {
  for(j in 1:N) {
    Sigma[i, j] = sigmaf^2 * exp(-(x[i] - x[j])^2 / (2 * l^2))
  }
}
    
(mvrnorm(n = 2, rep(0, N), Sigma) %>% t %>% as.tibble %>%
          mutate(x = x, V2 = exp(V2) / 5.0) %>%
          rename(mean = V1, sd = V2) %>%
    mutate(y = pmap(list(mean, sd), function(mean, sd) rnorm(10, mean, sd))) %>%
  (function(df) bind_cols(df, do.call("rbind", df$y) %>% as.tibble)) %>%
  select(-y) -> data) %>%
  select(-mean, -sd) %>%
  gather(variable, v, -x) %>%
  ggplot(aes(x, v)) +
  geom_point(aes(group = variable), alpha = 0.2) +
  geom_line(data = data %>%
            select(mean, sd, x) %>%
            gather(variable, v, -x), aes(colour = variable))

(fit = stan("models/heteroscedastic.stan", data = list(N = nrow(data),
                                                       y = data %>% select(V1:V5) %>% as.matrix,
                                                       x = data %>% pull(x),
                                                       M = 5,
                                                       sigmaf = sigmaf,
                                                       l = l), iter = 1000))

extract(fit, c("mu", "sigma")) %>% melt %>% dim

extract(fit, c("mu", "sigma")) %>% melt %>% as.tibble %>%
  rename(sample = iterations) %>%
  filter(sample < 50) %>%
  group_by(sample) %>% arrange(Var2) %>%
  spread(L1, value) %>% mutate(x = data$x) %>%
  select(-Var2) %>%
  gather(variable, v, c(-x, -sample)) %>%
  ggplot(aes(x, v)) +
  geom_line(aes(group = sample), alpha = 0.2) +
  geom_line(data = data %>%
              select(mean, sd, x) %>%
              rename(mu = mean,
                     sigma = sd) %>%
              gather(variable, v, -x), color = "red") +
  facet_grid(. ~ variable)
