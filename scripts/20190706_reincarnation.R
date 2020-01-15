#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
index <- as.numeric(Args[1])

# index <- 24
# first making a brain, from 0626_hpc.R -----------------------------------

source("./functions/functions_reports.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt


alphabeta <- c(1, 5, 0,
               0, 5, 1,
               0, 6, 0) %>% matrix(nrow = 3, byrow = TRUE) %>%
  as.data.frame()
colnames(alphabeta) <- c("p_low", "p_medium", "p_high")

l.alpha <- nrow(alphabeta)
nrounds <- 10
row_eps <-
  c(1:l.alpha) %>% rep(times = l.alpha) %>% rep(times = nrounds)
row_a <-
  c(1:l.alpha) %>% rep(each = l.alpha) %>% rep(times = nrounds)
r_ <- c(1:nrounds) %>% rep(each = l.alpha * l.alpha)
index_ <- 1:length(r_)

## to make the index file of the remaining rows
## (after adding (2,2) & (9,9)) uncomment the # %>%  ...
r_a_b <-
  data.frame(row_eps, row_a, r_, index_) # %>% filter(row_a>6|row_eps>6) %>% filter(round==1)

# r_a_b <- r_a_b %>% filter(row_eps==3 | row_a==3)

# making and saving the brain ---------------------------------------------

sampled.path <- "./data/5200-edges/"

brain_case <- partition_culture(
  round = r_a_b$r_[index],
  row_eps = r_a_b$row_eps[index],
  row_a = r_a_b$row_a[index],
  final.age = 10
)
save_vars(
  list.of.vars = "brain_case",
  prefix = paste(
    "life-01",
    brain_case@parameters$brain.code,
    brain_case@name,
    sep = "_"
  ),
  path = sampled.path
)

t0 <- Sys.time()
for (reincarnation in 2:10) {
  t1 <- Sys.time()
  
  pattern <- "_g-0.3k-5.2k"
  
  sampled.names <-
    list.files(path = sampled.path, pattern = "*.RData")
  
  # names <- sampled.names %>%
  #   my_gsub(paste0(".*",
  #                  pat.tmp,
  #                  "_"), "") %>%
  #   my_gsub("\\_.*", "") %>%
  #   sort()
  
  this.owner.oldest <-
    sampled.names[grepl(brain_case@name, sampled.names)] %>%
    sort() %>%
    tail(1)
  
  current.life <-
    this.owner.oldest %>% substring(6, 7) %>% as.numeric()
  
  # for(sampled in r.this){
  
  life.prefix <- sprintf("life-%02d", (current.life + 1))
  
  load(paste(sampled.path, this.owner.oldest, sep = ""))
  
  ## rebirth
  brain_case@initial <- brain_case@now
  # making initial$mat.con a list (standard form)
  brain_case@initial$mat.connectivity <-
    brain_case@now$mat.connectivity %>% list()
  # removing the long history of previous life
  brain_case@history <- brain_case@initial
  
  brain_case <- partition_culture(brain_case = brain_case,
                                  final.age = (brain_case@age$rewires %/% 1000) +
                                    100)
  
  Sys.time() - t1
  
  tryCatch({
    reports_netviz(brain_case)
  }, error = function(e) {
    print(paste("Error plotting", this.owner.oldest))
  })
  
  Sys.time() - t1
  
  
  save_vars(
    list.of.vars = "brain_case",
    prefix = paste(
      life.prefix,
      brain_case@parameters$brain.code,
      brain_case@name,
      sep = "_"
    ),
    path = sampled.path
  )
  
}

Sys.time() - t0
