library(MASS)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(deSolve)
library(GGally)
library(purrr)
library(parallel)
library(ggthemes)
library(tidyverse)
library(rstan)

setwd("~/gp")

N = 10

x = seq(-0.5, 0.5, length = N + 2)[2 : (N + 1)]
dx = x[[2]] - x[[1]]
xD = c(-0.5, x) + dx / 2.0

alpha = 1.0
l = 0.5

Sigma = matrix(0, nrow = N + 1, ncol = N + 1)
for(i in 1:(N + 1)) {
  for(j in 1:(N + 1)) {
    Sigma[i, j] = alpha^2 * exp(-(xD[i] - xD[j])^2 / (2 * l^2))
  }
}

inv_logit = function(x) { 1 / (1 + exp(-x)) }

D = exp(mvrnorm(1, rep(0, N + 1), Sigma))
plot(xD, D)

#D = rep(1.0, N + 1)#c(x, x[[N]] + dx) - dx / 2.0
#D[1] = 0.0
#D[N + 1] = 0.0

func = function(t, u, D) {
  dudt = rep(0, N)
  dudt[[1]] = (D[[2]] * (u[[2]] - u[[1]]) - D[[1]] * (u[[1]] - 1.0)) / dx^2
  dudt[[N]] = (D[[N + 1]] * (1.0 - u[[N]]) - D[[N]] * (u[[N]] - u[[N - 1]])) / dx^2
  for(i in 2 : (N - 1)) {
    dudt[[i]] = (D[[i + 1]] * (u[[i + 1]] - u[[i]]) - D[[i]] * (u[[i]] - u[[i - 1]])) / dx^2
  }
  list(dudt)
}

h = 0.01
times = seq(0, 0.1, by = h)

u0 = rep(0, N)
fout = inner_join(as_tibble(ode(y = u0, times = times, func = func, parms = D)[,]) %>%
                    gather(xi, u, 2:(N + 1)) %>%
                    mutate(xi = as.integer(xi)) %>%
                    mutate(unoise = u + rnorm(n(), sd = 0.05)),
  as_tibble(list(xi = 1:N, x = x)),
  by = "xi"
)

fout %>% filter(time == 0.1) %>%
  gather(name, u, c(u, unoise)) %>%
  ggplot(aes(x, u)) +
  geom_line(aes(group = name, color = name))

sdata = list(N = N,
             M = 5,
             t = 0.1,
             scale = 0.25,
             u0 = u0,
             u = (fout %>% filter(time == 0.1))$unoise,
             xD = xD)

fit = stan("models/diffusion_gp.stan", data = sdata, chains = 1, iter = 1000)

get_lines = function(fit, vnames, n = 100) {
  a = extract(fit, vnames)
  idxs = sample(nrow(a[[vnames[[1]]]]), n)
  
  out = as_tibble()
  for(i in 1:length(vnames)) {
    vname = vnames[[i]];
    d = a[[vname]][idxs,]
    colnames(d) <- 1:ncol(a[[vnames[[1]]]])
    d = as_tibble(d) %>% gather(time, data)
    d$name = vname
    
    out = bind_rows(out, d)
  }
  out %>% mutate(time = as.double(time))
}

out2 = inner_join(get_lines(fit, c("D"), 500) %>% rename(D = data) %>% select(time, D) %>% rename(xi = time),
                  as_tibble(list(xi = 1:(N + 1), xD = xD)),
                  by = "xi")

summary = out2 %>% group_by(xD) %>%
  summarize(mean = mean(D),
            m = quantile(D, 0.025),
            p = quantile(D, 0.975)) %>%
  ungroup()

summary %>% ggplot(aes(xD, mean)) +
  geom_ribbon(aes(ymin = m, ymax = p), alpha = 0.25, fill = "grey") +
  geom_line() +
  geom_line(aes(xD, p), alpha = 0.5) +
  geom_line(aes(xD, m), alpha = 0.5) +
  geom_line(data = as_tibble(list(xD = xD, D = D)), aes(xD, D), col = "red") +
  xlab("x") +
  ylab("Diffusion coefficient") +
  ggtitle("Diffusion coefficient (w/ est. 95% conf. intervals)")
