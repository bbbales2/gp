library(tidyverse)
library(lubridate)
library(rstan)

df0 = read_csv('~/Documents/nba/shots.csv', col_types = cols(.default = "?", GAME_CLOCK = "c"))

df = df0 %>% select(time = GAME_CLOCK, period = PERIOD, result = SHOT_RESULT, name = NAME) %>%
  filter((name == "Russell Westbrook" | name == "Kyrie Irving" | name == "Damian Lillard" | name == "Stephen Curry" | name == "John Wall") & period < 5) %>%
  mutate(resultBin = as.numeric(result == "made")) %>%
  mutate(time = ms(time)) %>%
  mutate(remaining_minutes = (4 - period) * 12 + as.numeric(as.duration(time) / 60.0))

df %>% ggplot(aes(remaining_minutes)) +
  geom_histogram(aes(y = ..density.., fill = factor(resultBin)), binwidth = 1.0, position = "fill") +
  geom_hline(yintercept = 0.5)# +
#facet_grid(period ~ .)

N = 500

df2 = df %>% mutate(gtime = 48 - remaining_minutes,
                    x = -(remaining_minutes - 24) / 48.0,
                    y = resultBin) %>%
  group_by(name) %>%
  sample_n(N) %>%
  ungroup() %>%
  arrange(name) %>%
  select(gtime, x, y, name)

names = (df2 %>% select(name) %>% unique())$name

sdata = list(N = N,
             M = 10,
             L = 5,
             scale = 0.25,
             x = df2$x,
             y = df2$y)

fit = stan('~/gp/models/multiple_players.stan', data = sdata, chains = 1)

get_lines = function(fit, vnames, n = 100) {
  a = extract(fit, vnames)
  idxs = sample(nrow(a[[vnames[[1]]]]), n)
  
  out = as_tibble()
  for(j in 1:length(names)) {
    for(i in 1:length(vnames)) {
      vname = vnames[[i]];
      d = a[[vname]][idxs,j,]
      colnames(d) <- (1:dim(a[[vnames[[1]]]])[[3]]) + (j - 1) * N
      d = as_tibble(d) %>% gather(itime, data)
      d$name = names[[j]]
      
      out = bind_rows(out, d)
    }
  }
  out %>% mutate(itime = as.double(itime))
}

a = get_lines(fit, c("f"), 500) %>% rename(f = data)
b = as_tibble(list(itime = 1:nrow(df2), time = df2$gtime))
out = inner_join(a, b, by = "itime")

summary = out %>% group_by(itime, name) %>%
  summarize(time = mean(time),
            mean = mean(f),
            m = mean - 2 * sd(f),
            p = mean + 2 * sd(f)) %>%
  ungroup()

summary %>% ggplot(aes(time, mean)) +
  #geom_ribbon(aes(ymin = m, ymax = p), alpha = 0.25) +
  geom_line(aes(group = name, color = name)) +
  geom_line(aes(time, p, group = name, color = name), alpha = 0.5, linetype = 'dashed') +
  geom_line(aes(time, m, group = name, color = name), alpha = 0.5, linetype = 'dashed') +
  #geom_point(data = out2, aes(time, f), size = 0.1, alpha = 0.01) +
  xlab("Game time") +
  ylab("Shooting percentage") +
  ggtitle("Hierarchical shooting percentage (w/ est. 95% conf. intervals)")

summary %>% arrange(time) %>% slice(c(1, n()))
