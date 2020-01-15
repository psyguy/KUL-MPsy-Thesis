rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt


# making the parameter file -----------------------------------------------

# levels_ <- c(0.5, 1, 1, 1, 2, 5)
# levels_ <- c(0.5, 1, 1, 2)
# levels.perm <- paste(rep(levels_, each = length(levels_)),
#                      rep(levels_, times = length(levels_)),
#                      sep = "x") %>% unique()
# 
# levels.perm2 <- paste(rep(levels.perm, each = length(levels.perm)),
#                       rep(levels.perm, times = length(levels.perm)),
#                       sep = "x") %>% unique()



# alpha <- c(0.5, 1, 1, 1, 2, 5, 2, 9)# %>% rep(times = nrounds)
# beta <- c(0.5, 1, 2, 5, 1, 1, 2, 9)# %>% rep(times = nrounds)
# alphabeta <- data.frame(alpha,beta)

alphabeta <- c(1,5,0,
               0,5,1,
               0,6,0) %>% matrix(nrow = 3, byrow = TRUE) %>%
  as.data.frame()
colnames(alphabeta) <- c("p_low", "p_medium", "p_high")

l.alpha <- nrow(alphabeta)
nrounds <- 3
row_eps <- c(1:l.alpha) %>% rep(times = l.alpha) %>% rep(times = nrounds)
row_a <- c(1:l.alpha) %>% rep(each = l.alpha) %>% rep(times = nrounds)
round <- c(1:nrounds) %>% rep(each = l.alpha*l.alpha) -1
index <- 1:length(round)

## to make the index file of the remaining rows
## (after adding (2,2) & (9,9)) uncomment the # %>%  ...
r_a_b <- data.frame(row_eps, row_a, round,index) # %>% filter(row_a>6|row_eps>6) %>% filter(round==1)


tt <- Sys.time()
for(r_a in 1:nrow(alphabeta)){
  for(r_eps in 1:nrow(alphabeta)){
    paste(r_a, r_eps) %>% print()
    # if(r_a<7 & r_eps<7) next
    brain_case <- partition_culture(round = 0,
                                    row_eps = r_eps,
                                    row_a = r_a,
                                    alphabeta = alphabeta,
                                    final.age = 1)
    save_vars(list.of.vars = "brain_case",
              prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))
  }
}
Sys.time() - tt




ab_pairs <- c("0.5x0.5")

round_alphabeta <- 0 %>% #c(1:6) %>%
  rep(each = length(levels.perm)) %>%
  paste(levels.perm, levels.perm, sep = "x") %>%
  as.data.frame()

# write.table(round_alphabeta,"round_alphabeta.txt",
#             row.names=FALSE,
#             col.names=FALSE)


# reading the round_alphabeta variable ------------------------------------

# r_eps_a <- read.table("round_alphabeta.txt")[1,] %>% as.character()

#!/usr/bin/env Rscript
# tt <- Sys.time()
# for(rab in round_alphabeta){
#   s <- round_alphabeta[rab,] %>% as.character()
#   t <- Sys.time()
#   s %>% print()
#   r_eps_a <-  (s %>% strsplit("x"))[[1]] %>% as.numeric()
#   
#   brain_case <- partition_culture(round = r_eps_a[1],
#                                   row_eps = r_eps_a[2],
#                                   row_a = r_eps_a[3],
#                                   final.age = 15)
#   save_vars(list.of.vars = "brain_case",
#             prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))
#   
#   paste(round_alphabeta[rab,], "took", (Sys.time()-t)) %>% print()
# }
# Sys.time() - tt

# making some small brains ------------------------------------------------

# paste(r_a, r_eps) %>% print()

