library(tidyverse)
library(rstan)

# simulate expt and x^2
expt_synth <- tibble(t = seq(-2,2,0.2)) %>% mutate(f = exp(t) + rnorm(nrow(.),sd = 0.05))
expt_synth %>% ggplot(aes(t,f)) + geom_point()

tsq_synth <- tibble(t = seq(-2,2,0.2)) %>% mutate(f = (1/2)*t^2 + rnorm(nrow(.),sd = 0.05))
tsq_synth %>% ggplot(aes(t,f)) + geom_point()


# estimate hyperparameters
fit_expt <- stan("stan/fit_hyperparameters.stan",
     data = list(N = nrow(expt_synth), t = expt_synth$t, y = expt_synth$f),
     chains = 1, iter = 1000)

fit_tsq <- stan("stan/fit_hyperparameters.stan",
                 data = list(N = nrow(tsq_synth), t = tsq_synth$t, y = tsq_synth$f),
                 chains = 1, iter = 1000)

get_ml_from_stan_samples <- function(fit) {
  s <- extract(fit)
  max_lik_idx <- which(s$lp__ == max(s$lp__))
  return(list(alpha = s$alpha[max_lik_idx],
              rho = s$rho[max_lik_idx],
              sigma = s$sigma[max_lik_idx]))
}

expt_hyperpar <- get_ml_from_stan_samples(fit_expt)
tsq_hyperpar <- get_ml_from_stan_samples(fit_tsq)

# get function GP estimates and derivative estimates
data_to_gp <- function(t, Xn, hyperpar, p, f) {
  
  # f is true function
  p_pars <- p(t = t, Xn = Xn, phi_n = hyperpar[1:2], sigma_n = hyperpar[[3]])
  gp <- tibble(t = t, mu = p_pars$mn, sigma = sqrt(diag(p_pars$Kn))) %>%
    mutate(f = f(t)) %>%
    mutate(lower = mu - 2*sigma, upper = mu + 2*sigma)
  
  return(gp)
}

Xn_posterior_expt <- data_to_gp(t = expt_synth$t, Xn = expt_synth$f, hyperpar = expt_hyperpar, p_Xn, exp)
Xn_posterior_expt %>% ggplot(aes(t, mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + geom_line(aes(y = f), color = "red")

dotXn_posterior_expt <- data_to_gp(t = expt_synth$t, Xn = expt_synth$f, hyperpar = expt_hyperpar, p_dotXn, exp)
dotXn_posterior_expt %>% ggplot(aes(t, mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + geom_line(aes(y = f), color = "red")

Xn_posterior_tsq <- data_to_gp(t = tsq_synth$t, Xn = tsq_synth$f, hyperpar = tsq_hyperpar, p_Xn, function(t) (1/2)*t^2)
Xn_posterior_tsq %>% ggplot(aes(t, mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + geom_line(aes(y = f), color = "red")

dotXn_posterior_tsq <- data_to_gp(t = tsq_synth$t, Xn = tsq_synth$f, hyperpar = tsq_hyperpar, p_dotXn, function(t) t)
dotXn_posterior_tsq %>% ggplot(aes(t, mu)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") + geom_line() + geom_line(aes(y = f), color = "red")
dotXn_posterior_tsq %>% mutate(ci = lower < f & f < upper) %>% summarize(ci = mean(ci))

#######################
#######################
# test create_p_dotXnS
#######################
#######################

# get mean and covariance of derivative GP at the time points we have data at
p <- p_dotXn(tn = expt_synth$t, Xn = expt_synth$f, phi_n = expt_hyperpar[1:2], sigma_n = expt_hyperpar[[3]])

# use the mean to estimate the hyper parameters of the GP regression
fit_dotXn <- stan("stan/fit_hyperparameters.stan",
                data = list(N = nrow(expt_synth), t = expt_synth$f, y = p$condMean),
                chains = 1, iter = 1000)

dotXn_expt_hyperpar <- get_ml_from_stan_samples(fit_dotXn)
qplot(expt_synth$f, p$condMean)

# get smoothed version of state as in the input in regression
p_state <- p_Xn(tn = expt_synth$t, Xn = expt_synth$f, phi_n = expt_hyperpar[1:2], sigma_n = expt_hyperpar[[3]])

# create sampling function
p_dotXnS <- create_p_dotXnS(Xn_list = list(p_state$condMean), mn = p$condMean, Kn = p$condVar, theta = dotXn_expt_hyperpar[1:2])

# sample from it
p_dotXnS(0.6)
p_dotXnS(1)
p_dotXnS(0.5)
p_dotXnS(0.1)
p_dotXnS(0.2)
p_dotXnS(1.2)

#
first_try <- function(xs) {
  p_dotXnS <- create_p_dotXnS(Xn_list = list(p_state$mn), mn = p$mn, Kn = p$Kn, theta = dotXn_expt_hyperpar[1:2])
  ret <- p_dotXnS(xs)
  tibble(mu = as.numeric(ret$mu), sigma = as.numeric(ret$sigma))
}

map(seq(0,4,0.1), first_try) %>% bind_rows %>% mutate(x = seq(0,4,0.1)) %>% mutate(lower = mu - 2*sigma, upper = mu + 2*sigma) %>%
  ggplot(aes(x,mu)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red")
