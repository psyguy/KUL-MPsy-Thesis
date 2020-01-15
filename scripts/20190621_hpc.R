rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt


# making the parameter file -----------------------------------------------

# levels_ <- c(0.5, 1, 1, 1, 2, 5)
# levels.perm <- paste(rep(levels_, each = length(levels_)),
#                rep(levels_, times = length(levels_)),
#                sep = "x") %>% unique()
# 
# round_alphabeta <- c(1:6) %>%
#   rep(each = length(levels.perm)) %>%
#   paste(levels.perm, sep = "x") %>%
#   as.data.frame()
# 
# write.table(round_alphabeta,"round_alphabeta.txt",
#             row.names=FALSE,
#             col.names=FALSE)


# reading the round_alphabeta variable ------------------------------------

# r_eps_a <- read.table("round_alphabeta.txt")[1,] %>% as.character()

#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
r_eps_a <- Args[1]

v1 <- gsub('"', '', r_eps_a)
r_eps_a <-  (v1 %>% strsplit("x"))[[1]]

# making some small brains ------------------------------------------------

# paste(r_a, r_eps) %>% print()
brain_case <- partition_culture(round = r_eps_a[1],
                                row_eps = r_eps_a[2],
                                row_a = r_eps_a[3],
                                final.age = 100)
save_vars(list.of.vars = "brain_case",
          prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))


