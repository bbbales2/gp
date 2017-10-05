library(MASS)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)

df_long_trim = readRDS("fake_data_proposal.rds")
(df_long_trim %>%
    filter(condition == "condition1", id == 1) %>%
    select(-condition, -id, -coordinate) %>%
    rename(x = time, y = position) %>%
    spread(trial, y) -> data) %>%
  gather(trial, y, -x) %>%
  ggplot(aes(x, y)) +
  geom_line(aes(group = trial), alpha = 0.2)

(fit = stan("models/heteroscedastic.stan", data = list(N = nrow(data),
                                                       y = data %>% select(-x) %>% as.matrix,
                                                       x = data %>% pull(x),
                                                       M = ncol(data) - 1), cores = 4, iter = 1000))

extract(fit, c("mu", "sigma")) %>% melt %>% dim

extract(fit, c("mu", "sigma")) %>% melt %>% as.tibble %>%
  rename(sample = iterations) %>%
  filter(sample < 50) %>% # only plot 50
  group_by(sample) %>% arrange(Var2) %>%
  spread(L1, value) %>% mutate(x = data$x) %>%
  select(-Var2) %>%
  gather(variable, v, c(-x, -sample)) %>%
  ggplot(aes(x, v)) +
  geom_line(aes(group = sample), alpha = 0.2) +
  facet_grid(. ~ variable)
