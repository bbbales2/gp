library(tidyverse)
library(rstan)

# simulate function
expt_synth <- tibble(t = seq(-2,2,0.04)) %>% mutate(f = exp(t) + rnorm(nrow(.),sd = 0.05))
expt_synth %>% ggplot(aes(t,f)) + geom_point()

# draw from a GP
gp_synth <- tibble(t = seq(-2,2,0.02))
K <- QQ(gp_synth$t, gp_synth$t, list(1,2))
gp_synth <- gp_synth %>% mutate(f = mvrnorm(1,mu = rep(0,nrow(gp_synth)),Sigma = K))# + rnorm(nrow(.),sd = 0.05))
gp_synth %>% ggplot(aes(t,f)) + geom_line()

# estimate hyperparameters
fit <- stan("stan/fit_hyperparameters.stan",
     data = list(N = nrow(expt_synth), t = expt_synth %>% pull(t), y = expt_synth %>% pull(f)),
     chains = 1, iter = 1000)

fit <- stan("stan/fit_hyperparameters.stan",
            data = list(N = nrow(gp_synth), t = gp_synth %>% pull(t), y = gp_synth %>% pull(f)),
            chains = 1, iter = 1000)

s <- extract(fit)
max_lik_idx <- which(s$lp__ == max(s$lp__))
s$alpha[max_lik_idx]
s$rho[max_lik_idx]
s$sigma[max_lik_idx]
qplot(s$alpha,s$rho) + geom_point(aes(x,y), color = "red", data = tibble(x=s$alpha[max_lik_idx],y=s$rho[max_lik_idx]))

# check if there's much uncertainty in the derivative that's attributable to the hyperparameters
f <- function(alpha,rho,sigma) {
  p <- p_dotXn(t = expt_synth %>% pull(t), Xn = expt_synth %>% pull(f), phi_n = list(alpha,rho), sigma_n = sigma)
  tibble(mn = p$mn %>% as.vector, sigma = diag(p$Kn))
}

dfs <- tibble(alpha = s$alpha, rho = s$rho, sigma = s$sigma) %>% pmap(f) %>%
  map2(1:500, function(x,y) x %>% mutate(lower = mn - 2*sigma, upper = mn + 2*sigma, sample = y, t = expt_synth %>% pull(t))) %>% bind_rows
dfs %>% ggplot(aes(t,mn,group = sample)) + geom_line(alpha = 0.1) + stat_function(fun = exp, color = "red")

dfs %>% group_by(t) %>% summarize(med = median(mn), lower = min(lower), upper = max(upper)) %>%
  ggplot(aes(t,med)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line() +
  stat_function(fun = function(x) x, color = "red")

# get mn and Kn now that we have hyperparameter point estimates
p_Xn_pars <- p_Xn(t = expt_synth$t, Xn = expt_synth$f, phi_n = list(alpha = 3.24, rho = 1.29), sigma_n = 0.055)
mn <- p_Xn_pars$mn %>% as.vector
Kn <- p_Xn_pars$Kn
Xn_est <- tibble(t = expt_synth$t, mu = mn, sigma = sqrt(diag(Kn))) %>%
  mutate(f = exp(t)) %>% mutate(lower = mu - 2*sigma, upper = mu + 2*sigma)# %>%
  #ggplot(aes(t,mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + stat_function(fun = exp, color = "red")

p_dotXn_pars <- p_dotXn(t = expt_synth %>% pull(t), Xn = expt_synth %>% pull(f), phi_n = list(alpha = 3.24, rho = 1.29), sigma_n = 0.055)
mn <- p_dotXn_pars$mn %>% as.vector
Kn <- p_dotXn_pars$Kn

dotXn_est <- tibble(t = expt_synth %>% pull(t), mu = mn, sigma = diag(Kn)) %>% mutate(f = exp(t)) %>% mutate(lower = mu - 2*sqrt(sigma), upper = mu + 2*sqrt(sigma))
dotXn_est %>% ggplot(aes(t,mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + stat_function(fun = exp, color = "red")
