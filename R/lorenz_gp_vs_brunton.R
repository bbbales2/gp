library(tidyverse)
library(deSolve)

#define ode for deSolve
Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dx <- sigma*(y-x)
    dy <- x*(rho-z) - y
    dz <- x*y - beta*z
    
    # return the rate of change
    list(c(dx, dy, dz))
  })
}

#use same parameters as Brunton et al.
sim_lorenz <- function(h, t = 20, theta = c(sigma = 10, beta = 8/3, rho = 28), y0 = c(x = -8, y = 7, z = 27)) {
  ode(y = y0, times = seq(0, t, by = h), func = Lorenz, parms = theta) %>%
    as.data.frame %>%
    as_tibble
}

#simulate and plot
lorenz_data <- sim_lorenz(h = 0.1)
lorenz_data %>%
  gather(state,value,-time) %>%
  ggplot(aes(time,value)) +
  geom_point() +
  geom_line() +
  facet_grid(state ~ .)

#
