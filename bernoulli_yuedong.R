library(MASS)
library(car)
library(assist)
library(tidyverse)
library(rstan)
library(ggplot2)
library(shinystan)

setwd("~/gp")

N = 200
x = seq(0.0, 1.0, length = N)
y = 10.0 * cos(x * 2.7) + 1.0 + rnorm(N) * 0.5;

qplot(x, y)

f3 = ssr(y ~ x, rk = cubic(x), scale = cbind(1, x))
as_tibble(list(x = x, y = y, mu = predict(f3)$fit, std = predict(f3)$pstd)) %>%
  mutate(p = mu + 2.0 * std, m = mu - 2.0 * std) %>%
  ggplot(aes(x, mu)) +
  geom_line() +
  geom_line(aes(x, m), linetype = "dashed") +
  geom_line(aes(x, p), linetype = "dashed") +
  geom_point(aes(x, y))

fit = stan("models/bernoulli_yuedong2.stan", data = list(N = N,
                                                        x = x,
                                                        y = y), chains = 1, cores = 4, iter = 1000)

sm = stan_model("models/bernoulli_yuedong.stan")

fit2 = optimizing(sm, data = list(N = N,
                                  x = x,
                                  y = y))

a = as_tibble(extract(fit, "yhat")$yhat)
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
  geom_point(data = as_tibble(list(x = x, y = y)), aes(x, y), colour = "purple") +
  xlab("x") +
  ylab("y")

{
  k0 = function(x) {
    1.0
  }
  
  k1 = function(x) {
    x - 0.5
  }
  
  k2 = function(x) {
    k1_ = k1(x)
    0.5 * (k1_ * k1_ - 1.0 / 12.0)
  }
  
  k4 = function(x) {
    k1_ = k1(x)
    (1.0 / 24.0) * (k1_ * k1_ * k1_  * k1_ - 0.5 * k1_ * k1_ + 7.0 / 240.0)
  }
}

Sigma = matrix(0, nrow = N, ncol = N);
vk0 = k0(x)
vk1 = k1(x)

for(i in 1:N) {
  for(j in 1:N) {
    Sigma[i, j] = k2(x[i]) * k2(x[j]) - k4(abs(x[i] - x[j]));
  }
}

for(i in 1:N) {
  Sigma[i, i] = Sigma[i, i];
  
  vk0[i] = k0(x[i]);
  vk1[i] = k1(x[i]);
}


qplot(x, mvrnorm(mu = rep(0.0, N), Sigma = Sigma))

y = k0(x) + k2(x)

qplot(x, y)

a = vk0 %*% y / N
b = vk1 %*% y / N
cc = t(y) %*% Sigma

plot(x, cc + a * vk0 + b * vk1 - y)
plot(x, cc)
L = t(chol(Sigma));

L = cholesky_decompose(Sigma);