library(tidyverse)
library(lubridate)
library(rstan)
library(shinystan)

#df2 = read_csv('~/gp/discourse_westbrook/westbrook.csv')

df0 = read_csv('~/Documents/nba/shots.csv', col_types = cols(.default = "?", GAME_CLOCK = "c"))

df = df0 %>% select(time = GAME_CLOCK, period = PERIOD, result = SHOT_RESULT, name = NAME, pts_type = PTS_TYPE) %>%
  filter(name == "Russell Westbrook" & period < 5) %>%
  mutate(resultBin = as.numeric(result == "made")) %>%
  mutate(time = ms(time)) %>%
  mutate(remaining_minutes = (4 - period) * 12 + as.numeric(as.duration(time) / 60.0))

df2 = df %>% mutate(gtime = 48 - remaining_minutes,
                    x = -(remaining_minutes - 24) / 48.0,
                    y = resultBin) %>%
  select(gtime, period, x, y, result) %>% sample_n(100)

sdata = list(N = dim(df2)[[1]],
             M = 10,
             scale = 0.25,
             x = df2$x,
             y = df2$y)

N = dim(df2)[[1]]
M = 10
scale = 0.25
x = df2$x
y = df2$y

stan_rdump(c('N', 'M', 'scale', 'x', 'y'), file = '~/gp/sdata.dat')

fit = stan('~/gp/models/westbrook.stan', data = sdata, chains = 1, cores = 1, iter = 2000)

launch_shinystan(fit)

sdata = list(N = dim(df2)[[1]],
             T = 4,
             x = df2$period,
             y = df2$y)

fit2 = stan('~/gp/models/westbrook_grouped.stan', data = sdata, chains = 4, cores = 4)

get_lines = function(fit, vnames, n = 100) {
  a = extract(fit, vnames)
  idxs = sample(nrow(a[[vnames[[1]]]]), n)
  
  out = as_tibble()
  for(i in 1:length(vnames)) {
    vname = vnames[[i]];
    d = a[[vname]][idxs,]
    colnames(d) <- 1:ncol(a[[vnames[[1]]]])
    d = as_tibble(d) %>% gather(idx, data)
    d$name = vname
    
    out = bind_rows(out, d)
  }
  out %>% mutate(idx = as.numeric(idx))
}

a = get_lines(fit, c("f"), 50) %>% rename(f = data) %>% select(idx, f)
b = as_tibble(list(idx = 1:nrow(df2), time = df2$gtime))
out2 = inner_join(a, b, by = "idx")

a = get_lines(fit, c("f"), 1000) %>% rename(f = data) %>% select(idx, f)
b = as_tibble(list(idx = 1:nrow(df2), time = df2$gtime))
out = inner_join(a, b, by = "idx")

summary = out %>% group_by(time) %>%
          summarize(mean = mean(f),
                    m = quantile(f, 0.025),
                    p = quantile(f, 0.975),
                    sd = sd(f),
                    esd = sqrt(mean * (1 - mean))) %>%
          ungroup()

out3 = extract(fit2, "f")$f
colnames(out3) = 1:4
summary2 = left_join(as_tibble(out3) %>% gather(period, f),
          as_tibble(list(period = c("1", "2", "3", "4"), time = c(6.0, 18.0, 30.0, 42.0))),
          by = "period") %>% 
  group_by(time) %>%
  summarize(mean = mean(f),
            m = quantile(f, 0.025),
            p = quantile(f, 0.975)) %>%
  ungroup()


summary %>% ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = m, ymax = p), alpha = 0.25, fill = "blue") +
  geom_line() +
  geom_line(aes(time, p), col = "blue", alpha = 0.5) +
  geom_line(aes(time, m), col = "blue", alpha = 0.5) +
  geom_point(data = out2, aes(time, f), size = 0.1, alpha = 0.01) +
  #geom_boxplot(data = summary2, aes(time, f, group = time)) +
  geom_errorbar(data = summary2, aes(time, ymin = m, ymax = p)) +
  geom_point(data = summary2, aes(time, mean)) +
  xlab("Game time") +
  ylab("Shooting percentage") +
  ggtitle("Russell Westbrook's shooting percentage (w/ est. 95% conf. intervals)")

summary %>% arrange(time) %>% slice(c(1, n()))
