rm(list = ls())

source("./functions/functions_trial.R")


# time series of nodes ----------------------------------------------------

library(ggplot2)
library(reshape)


# data <- data.frame(time = seq(0, 23), noob = rnorm(24), plus = runif(24), extra = rpois(24, lambda = 1))
# Molten <- melt(data, id.vars = "time")
# ggplot(Molten, aes(x = time, y = value, colour = variable)) + geom_line()



tmp_seq <- seq(1, 9001, 10)

sbs_h <- (brain_case@history$activities[tmp_seq,1:30]) %>% as.data.frame()

# sbs_h[,12] %>% plot()

data <- sbs_h %>% cbind(tmp_seq)

# Molten$variable %>% str() 

# data <- data.frame(time = seq(0, 23), noob = rnorm(24), plus = runif(24), extra = rpois(24, lambda = 1))
Molten <- melt(data, id.vars = "tmp_seq")
ggplot(Molten, aes(x = tmp_seq, y = value, colour = variable)) + geom_line()


# network measures --------------------------------------------------------

g <- graph.ring(10)
shortest.paths(g)
get.shortest.paths(g, 5)
average.path.length(g)
