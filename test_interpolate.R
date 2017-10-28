library(tidyverse)
library(ggplot2)
library(rstan)
library(Rcpp)
library(shinystan)

sourceCpp("covariance.cpp")

cpp_src = stanc("models/cubic_interpolated_gp.stan", allow_undefined = TRUE)$cppcode
cpp_src_split = strsplit(cpp_src, "\n")[[1]]
first_match = grep("approx_Lz", cpp_src_split)[[1]]
cat(cpp_src_split[(first_match - 2) : (first_match + 4)], sep = "\n")

model1 = stan_model("models/cubic_interpolated_gp.stan",
                   allow_undefined = TRUE,
                   includes = paste0('\n#include "',
                                     file.path(getwd(),
                                               'models/cubic_interpolated_gp.hpp'), '"\n'))

model2 = stan_model("models/exact_gp.stan")

i = 2
N = i * 10
x = seq(0.0, 10.0, length = N)
y = sin(x) + rnorm(N, 0.0, 0.2)
P = 3
lp = seq(qgamma(0.01, 4.0, 4.0), qgamma(0.99, 4.0, 4.0), length = P)

Ls = list()
dLdls = list()
for(p in 1:P) {
  out = rbf_cov_chol(x, lp[p])
  Ls[[p]] = out$L
  dLdls[[p]] = out$dLdl
}

#list(x = x, y = y) %>%
#  as.tibble %>%
#  ggplot(aes(x, y)) +
#  geom_point()

fit1 = sampling(model1,
                data = list(N = N, x = x, y = y, P = P, lp = lp, Ls = Ls, dLdls = dLdls),
                iter = 1000, chains = 1)

fit2 = sampling(model2,
                data = list(N = N, x = x, y = y),
                iter = 1000, chains = 1)

bind_rows(extract(fit1, c("l", "sigma")) %>% as.tibble %>% mutate(which = "approx"),
          extract(fit2, c("l", "sigma")) %>% as.tibble %>% mutate(which = "exact")) %>%
  gather(parameter, value, c("l", "sigma")) %>%
  ggplot(aes(value)) +
  geom_density(aes(fill = which), alpha = 0.5) +
  facet_grid(. ~ parameter)

bind_rows(extract(fit1, c("f"))$f %>% as.tibble %>% setNames(x) %>%
            gather(x, value) %>% mutate(x = as.numeric(x)) %>% mutate(which = "approx"),
          extract(fit2, c("f"))$f %>% as.tibble %>% setNames(x) %>%
            gather(x, value) %>% mutate(x = as.numeric(x)) %>% mutate(which = "exact")) %>%
  group_by(x, which) %>%
  summarize(q1 = quantile(value, 0.025),
            q2 = quantile(value, 0.25),
            q3 = quantile(value, 0.75),
            q4 = quantile(value, 0.975)) %>%
  ungroup() %>%
  ggplot(aes(x)) +
  geom_ribbon(aes(ymin = q2, ymax = q3, fill = which), alpha = 0.5) +
  geom_line(aes(y = q1, colour = which), linetype="dashed") +
  geom_line(aes(y = q4, colour = which), linetype="dashed") +
  ylab("Gaussian process output") +
  ggtitle("50% intervals solid,\n95% intervals dashed.\napprox vs. exact")
  

timings = tibble(N = numeric(), t = numeric(), which = character())
for(i in 1:5) {
  for(j in 1:3) {
    N = i * 10
    x = seq(0.0, 10.0, length = N)
    y = sin(x) + rnorm(N, 0.0, 0.2)
    P = 10
    lp = seq(qgamma(0.01, 4.0, 4.0), qgamma(0.99, 4.0, 4.0), length = 10)
    
    Ls = list()
    dLdls = list()
    for(p in 1:P) {
      out = rbf_cov_chol(x, lp[p])
      Ls[[p]] = out$L
      dLdls[[p]] = out$dLdl
    }
    
    #list(x = x, y = y) %>%
    #  as.tibble %>%
    #  ggplot(aes(x, y)) +
    #  geom_point()
    
    t1 = as.numeric(system.time({fit1 = sampling(model1,
                                data = list(N = N,
                                            x = x,
                                            y = y,
                                            P = P,
                                            lp = lp,
                                            Ls = Ls,
                                            dLdls = dLdls),
                     iter = 1000,
                     chains = 1)}))[1]
    
    timings = timings %>% add_row(N = N, t = t1, which = 'approx')
    
    t2 = as.numeric(system.time({fit2 = sampling(model2,
                 data = list(N = N,
                             x = x,
                             y = y),
                 iter = 1000,
                 chains = 1)}))[1]
    
    timings = timings %>% add_row(N = N, t = t2, which = 'exact')
  }
}

timings %>%
  ggplot(aes(N, t)) +
  geom_point(aes(color = which)) +
  ylab("Runtime (s)") +
  ggtitle("Approx vs. Exact runtime for 500 \nwarmup 500 samples 1 chain")
