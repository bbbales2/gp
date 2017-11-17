library(tidyverse)
library(ggplot2)

y1 = 5.0
y2 = 2.0
k1 = -5.0
k2 = 3.0
x1 = 1.0
x2 = 1.75

x = seq(x1 - 0.1, x2 + 0.1, length = 100)

t = (x - x1) / (x2 - x1)

a = k1 * (x2 - x1) - (y2 - y1)
b = -k2 * (x2 - x1) + (y2 - y1)

q = (1 - t) * y1 + t * y2 + t * (1 - t) * (a * (1 - t) + b * t)

qplot(x, q)
