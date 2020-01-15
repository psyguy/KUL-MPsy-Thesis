#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
index <- as.numeric(Args[1])


#rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt


alphabeta <- c(1,5,0,
               0,5,1,
               0,6,0) %>% matrix(nrow = 3, byrow = TRUE) %>%
  as.data.frame()
colnames(alphabeta) <- c("p_low", "p_medium", "p_high")

l.alpha <- nrow(alphabeta)
nrounds <- 10
row_eps <- c(1:l.alpha) %>% rep(times = l.alpha) %>% rep(times = nrounds)
row_a <- c(1:l.alpha) %>% rep(each = l.alpha) %>% rep(times = nrounds)
r_ <- c(1:nrounds) %>% rep(each = l.alpha*l.alpha)
index_ <- 1:length(r_)

## to make the index file of the remaining rows
## (after adding (2,2) & (9,9)) uncomment the # %>%  ...
r_a_b <- data.frame(row_eps, row_a, r_, index_) # %>% filter(row_a>6|row_eps>6) %>% filter(round==1)

# r_a_b <- r_a_b %>% filter(row_eps==3 | row_a==3)

# making and saving the brain ---------------------------------------------

brain_case <- partition_culture(round = r_a_b$r_[index],
                                row_eps = r_a_b$row_eps[index],
                                row_a = r_a_b$row_a[index],
                                final.age = 100)
save_vars(list.of.vars = "brain_case",
          prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))


