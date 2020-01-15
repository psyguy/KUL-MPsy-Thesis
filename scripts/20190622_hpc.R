#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
index <- as.numeric(Args[1])


#rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt


alpha <- c(0.5, 1, 1, 1, 2, 5)# %>% rep(times = nrounds)
beta <- c(0.5, 1, 2, 5, 1, 1)# %>% rep(times = nrounds)
alphabeta <- data.frame(alpha,beta)

l.alpha <- nrow(alphabeta)
nrounds <- 3
row_eps <- c(1:l.alpha) %>% rep(times = l.alpha) %>% rep(times = nrounds)
row_a <- c(1:l.alpha) %>% rep(each = l.alpha) %>% rep(times = nrounds)
round <- c(1:nrounds) %>% rep(each = l.alpha*l.alpha)

r_a_b <- data.frame(row_eps, row_a, round)


# making and saving the brain ---------------------------------------------

brain_case <- partition_culture(round = r_a_b$round[index],
                                row_eps = r_a_b$row_eps[index],
                                row_a = r_a_b$row_a[index],
                                final.age = 100)
save_vars(list.of.vars = "brain_case",
          prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))


