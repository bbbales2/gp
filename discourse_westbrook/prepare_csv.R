library(tidyverse)
library(lubridate)
library(rstan)
library(shinystan)

df0 = read_csv('~/Documents/nba/shots.csv', col_types = cols(.default = "?", GAME_CLOCK = "c"))

df = df0 %>% select(time = GAME_CLOCK, period = PERIOD, result = SHOT_RESULT, name = NAME, pts_type = PTS_TYPE) %>%
  filter(name == "Russell Westbrook" & period < 5) %>%
  mutate(resultBin = as.numeric(result == "made")) %>%
  mutate(time = ms(time)) %>%
  mutate(remaining_minutes = (4 - period) * 12 + as.numeric(as.duration(time) / 60.0))

df2 = df %>% mutate(gtime = 48 - remaining_minutes,
                    x = -(remaining_minutes - 24) / 48.0,
                    y = resultBin) %>%
  select(gtime, period, x, y, result)

write_csv(df2, '~/gp/discourse_westbrook/westbrook.csv')
